"""
MCFLIRT-style framewise head-motion estimation.

Each volume of a 4D series is rigidly registered (6 dof, Euler3DTransform) to a
fixed 3D reference, following the logic of FSL MCFLIRT's default estimation:
  * multi-stage schedule: each stage resamples the *reference only* to an
    isotropic grid (moving volumes stay native and unsmoothed) and runs one
    sweep of 1D line searches over rx, ry, rz, tx, ty, tz, each stopping when
    its bracket is narrower than one tolerance
  * each stage has its own rotation (rad) and translation (mm) tolerances
  * stage 1 is chained: the middle volume is registered from identity, then
    each volume starts from its neighbour's result, outward in both directions;
    later stages start each volume from its own previous-stage estimate
  * schedule sizes in mm (grid spacings, translation tolerances) are scaled 
    according to brain size of the input data (estimated from an input mask 
    file), scaling is assessed relative to human brain size (which is the 
    reference scale for default parameters)

Differences from MCFLIRT that remain by design: the chain starts at the middle
volume (MCFLIRT with -reffile starts at volume 0); the stage reference is
resampled with BSpline5 and rounded grid sizes (FSL: trilinear point sampling,
truncated sizes); the rotation centre is the reference's centre of mass; FSL's
fuzzy FOV-edge weighting in the cost is not reproduced. The sitk
backend also uses ITK's Brent / golden-ratio line search and ITK's correlation
metric instead of FSL's.
"""
import os
import numbers
import concurrent.futures
import functools

import numpy as np
import SimpleITK as sitk
from tqdm import tqdm

# 151mm spherical diameter estimated from the MNI's brain mask,
# provides a reference for the default scale of MCFLIRT's parameters
HUMAN_BRAIN_SIZE_REF = 151

# 9.1mm spherical diameter estimated from the DSURQE template mask
MOUSE_BRAIN_SIZE_REF = 9.10

MOUSE_TO_HUMAN_SCALING = HUMAN_BRAIN_SIZE_REF / MOUSE_BRAIN_SIZE_REF  # ~16.6X

# Length factor that makes a human-scale schedule, once scaled to a mouse brain,
# equal to MCFLIRT run on x10 header-scaled mouse data (~1.66)
_X10 = MOUSE_TO_HUMAN_SCALING / 10

# Per stage: (grid spacing in mm, rotation tolerance in rad, translation
# tolerance in mm), at human brain scale. Spacings and translation tolerances
# are multiplied by brain_scale at run time; rotation tolerances are not.
# MCFLIRT's tolerances are set_param_tols (0.005 rad, 0.2 mm) times the stage
# factor (0.8, 0.8, 0.1).
SCHEDULES = {
    # MCFLIRT's hard-coded default 3-stage schedule
    'fine': [
        (8.0, 0.004, 0.16),
        (4.0, 0.004, 0.16),
        (4.0, 0.0005, 0.02),
    ],
    # MCFLIRT with lengths x _X10: for a mouse brain, spacings and translation
    # tolerances match MCFLIRT run on x10 header-scaled data (0.8 / 0.4 / 0.4 mm);
    # rotation tolerances are MCFLIRT's, as header scaling leaves angles unchanged
    'coarse': [
        (8.0 * _X10, 0.004, 0.16 * _X10),
        (4.0 * _X10, 0.004, 0.16 * _X10),
        (4.0 * _X10, 0.0005, 0.02 * _X10),
    ],
    # the first stage can propagate errors to neighboring slices because of the
    # transform chain. A finer first stage minimize this risk.
    # The later stages can retain 'coarse' for efficiency, and avoiding overfitting
    # that can arise at finer scales.
    'balanced': [
        (8.0, 0.004, 0.16),
        (4.0 * _X10, 0.004, 0.16), # the translation is not scaled so that we do not increase tolerance at later stages
        (4.0 * _X10, 0.0005, 0.02 * _X10),
    ],
    # 'balanced' plus a finer 4th stage
    'stringent': [
        (8.0, 0.004, 0.16),
        (4.0 * _X10, 0.004, 0.16),
        (4.0 * _X10, 0.0005, 0.02 * _X10),
        (2.0 * _X10, 0.00035, 0.014 * _X10),
    ],
}


def _as_float_image(img):
    if isinstance(img, sitk.Image):
        if img.GetPixelID() == sitk.sitkFloat32:
            return img
        if img.GetDimension() <= 3:
            return sitk.Cast(img, sitk.sitkFloat32)
        # sitk.Cast does not support 4D images: go through numpy
        arr = sitk.GetArrayFromImage(img).astype(np.float32)
        out = sitk.GetImageFromArray(arr, isVector=False)
        out.CopyInformation(img)
        return out
    elif os.path.isfile(img):
        return sitk.ReadImage(img, sitk.sitkFloat32)  # works for 4D
    raise ValueError(f"{img} is neither a file nor an SITK image.")


def _read_schedule(schedule):
    if isinstance(schedule, str):
        if schedule not in SCHEDULES:
            raise ValueError(f"{schedule} is not among the pre-built schedules: {list(SCHEDULES.keys())}")
        return SCHEDULES[schedule]
    error = ValueError(f"{schedule} doesn't have the expected schedule format: "
                       "a non-empty list of (spacing mm, rotation tol rad, translation tol mm) "
                       "tuples of positive numbers.")
    if not isinstance(schedule, (list, tuple)) or len(schedule) == 0:
        raise error
    for stage in schedule:
        if not isinstance(stage, (list, tuple)) or len(stage) != 3:
            raise error
        for param in stage:
            if not isinstance(param, numbers.Real) or isinstance(param, bool) or param <= 0:
                raise error
    return [tuple(stage) for stage in schedule]


def brain_size_mm(brain_mask_file):
    # given the brain mask volume, estimate the brain spherical diameter
    mask = sitk.ReadImage(brain_mask_file, sitk.sitkUInt8) > 0   # binary, label 1
    stats = sitk.LabelShapeStatisticsImageFilter()
    stats.Execute(mask)
    return 2 * stats.GetEquivalentSphericalRadius(1)


def isotropic_upsample_and_pad(
    image, interpolation=sitk.sitkBSpline5, clip_negative=True, spacing=None, pad=2
):
    """
    Resample the image to isotropic spacing.

    Parameters:
        image (sitk.Image): Input SimpleITK image.
        interpolation: SimpleITK interpolation method
                       (e.g., sitk.sitkLinear, sitk.sitkNearestNeighbor, sitk.sitkBSpline).
        clip_negative (bool): Set negative values (interpolation overshoot) to 0.
        spacing (float or None): Target isotropic spacing. If None, uses the
                       smallest existing spacing (original behaviour).
        pad (int): Number of padding, i.e. duplicate of the outer slice on every side.

    Returns:
        sitk.Image: Isotropically resampled image.
    """
    original_spacing = image.GetSpacing()
    if spacing is None:
        spacing = min(original_spacing)
    if np.allclose(original_spacing, image.GetDimension() * (spacing,)):
        # the image is already at the target isotropic spacing
        resampled_image = image
    else:
        # Compute new size to maintain physical extent
        original_size = image.GetSize()
        new_size = [
            int(round(original_size[i] * original_spacing[i] / spacing))
            for i in range(len(original_spacing))
        ]

        # Resample
        resampled_image = sitk.Resample(
            image,
            new_size,
            sitk.Transform(),  # Identity transform
            interpolation,
            image.GetOrigin(),
            (spacing,) * len(original_spacing),
            image.GetDirection(),
            0,  # Default pixel value
            image.GetPixelID(),
        )
        if clip_negative:
            resampled_image = resampled_image * sitk.Cast(
                resampled_image > 0, resampled_image.GetPixelID()
            )
    if pad > 0:
        # Duplicate the outer slice of the image `pad` times on every side.
        dim = image.GetDimension()
        for i in range(pad):
            resampled_image = sitk.MirrorPad(resampled_image, [1] * dim, [1] * dim)
    return resampled_image


def register_pair_sitk_mcflirt(
    fixed,
    moving,
    initial_transform,
    tol,
    metric="correlation",
    sweeps=1,
    boundguess=10.0,
    num_threads=1,
):
    """One MCFLIRT-style stage: Powell line searches in units of `tol`.

    fixed             : reference already resampled to the stage grid
    moving            : native 3D volume
    initial_transform : Euler3DTransform (its centre is kept)
    tol               : 6 tolerances (rx, ry, rz in rad; tx, ty, tz in mm)
    metric            : "correlation" (MCFLIRT's normcorr) or "mattes"
    sweeps            : 1 = MCFLIRT; >1 enables full Powell iterations
    boundguess        : initial line-search bracket, in tolerances

    ITK mapping: numberOfIterations = sweeps - 1 (ITK's loop is `<=`), and a
    huge valueTolerance for a single sweep so ITK returns before Powell's extra
    line search. Line units are tol / c, so stepLength = boundguess * c starts
    the bracket at `boundguess` tolerances and stepTolerance = c / 2 stops it
    at one tolerance; c is small so ITK's relative Brent tolerance
    (stepTolerance * |x|) stays negligible next to the absolute one.
    """
    R = sitk.ImageRegistrationMethod()
    if metric == "correlation":
        R.SetMetricAsCorrelation()
    elif metric == "mattes":
        R.SetMetricAsMattesMutualInformation(numberOfHistogramBins=32)
    else:
        raise ValueError(f"unknown metric {metric}")
    R.SetMetricSamplingStrategy(R.NONE)
    R.SetInterpolator(sitk.sitkLinear)
    R.SetShrinkFactorsPerLevel([1])
    R.SetSmoothingSigmasPerLevel([0.0])

    c = 0.02  # line units = tol / c
    R.SetOptimizerAsPowell(
        numberOfIterations=sweeps - 1,
        maximumLineIterations=100,
        stepLength=boundguess * c,
        stepTolerance=c / 2,
        valueTolerance=1e10 if sweeps == 1 else 1e-8,
    )
    R.SetOptimizerScales([c / t for t in tol])
    R.SetInitialTransform(initial_transform, inPlace=False)
    R.SetNumberOfThreads(num_threads)

    try:
        return R.Execute(fixed, moving).GetBackTransform()
    except RuntimeError as e:
        print(f"Registration failed ({e}); keeping the initial transform")
        return sitk.Euler3DTransform(initial_transform)


def register_pair_fsl_mcflirt(
    fixed,
    moving,
    initial_transform,
    tol,
    metric="correlation",
    sweeps=1,
    boundguess=10.0,
):
    """Same inputs/outputs as register_pair_sitk_mcflirt, using FSL's optimiser.

    Each sweep runs a 1D line search along rx, ry, rz, tx, ty, tz in turn.
    The first sweep brackets at `boundguess` tolerances; later sweeps use 1
    tolerance (MCFLIRT's boundguess = [10, 1]) and stop early once no parameter
    moved by more than its tolerance. Only metric="correlation" is implemented.
    """
    if metric != "correlation":
        raise ValueError("register_pair_fsl_mcflirt only implements metric='correlation'")
    center = initial_transform.GetCenter()
    cost = _normcorr_cost(fixed, moving, center)
    p = np.array(initial_transform.GetParameters(), dtype=float)
    fval = cost(p)
    for sweep in range(sweeps):
        bg = boundguess if sweep == 0 else 1.0
        start = p.copy()
        for n in range(6):
            def f1d(x, n=n):
                q = p.copy()
                q[n] += x
                return cost(q)

            x, fval = _optimise1d(f1d, tol[n], fval, bg)
            p[n] += x
        if sweep > 0 and np.all(np.abs(p - start) < tol):
            break

    out = sitk.Euler3DTransform()
    out.SetCenter(center)
    out.SetParameters(tuple(float(v) for v in p))
    return out


def framewise_register_mcflirt_like(
    moving_img,
    ref_img,
    schedule='balanced',
    brain_scale=None,
    brain_size_mask=None,
    metric="correlation",
    padding=0,
    sweeps=1,
    parallel=False,
    max_workers=os.cpu_count(),
    backend="fsl_mcflirt",
    verbose=True,
    only_print_schedule=False,
):
    """Estimate rigid head motion for every volume of a 4D series against a
    fixed 3D reference, with MCFLIRT's multi-stage strategy.

    Algorithm
    ---------
    For each stage (spacing, rotation tol, translation tol) of the schedule:
      1. The reference is resampled to an isotropic grid of `spacing * brain_scale`mm. 
      Moving volumes are used at native resolution, without smoothing.
      2. Each volume is registered to that grid with one optimiser call from
         `backend`, with the stage's rotation tolerance on rx, ry, rz and its
         translation tolerance, multiplied by `brain_scale`, on tx, ty, tz.
      3. Initialisation: in stage 1 the middle volume starts from identity and
         the others are chained outward from it (each starts from its
         neighbour's stage-1 result, in both directions). In later stages each
         volume starts from its own previous-stage estimate.

    Parameters
    ----------
    moving_img : sitk.Image or str
        4D series (image or file path). Cast to float32.
    ref_img : sitk.Image or str
        3D registration target (image or file path)
    schedule : str or list of (float, float, float)
        Name of a pre-built schedule in SCHEDULES ('coarse', 'fine', 'balanced', 
        'stringent') or a list of (grid spacing in mm, rotation tolerance in
        rad, translation tolerance in mm) tuples, at human brain scale.
    brain_scale : float or None
        Brain size / human brain size. Grid spacings and
        translation tolerances are multiplied by it; rotation tolerances are
        not. If None, estimated from `brain_size_mask`, else 1 (human scale).
    brain_size_mask : str or None
        Brain mask file (e.g. of the reference template) whose equivalent
        spherical diameter gives `brain_scale`. Must be in true mm.
    metric : str
        "correlation" or "mattes"
    padding : int
        Number of duplicated edge slices added around each stage reference.
    sweeps : int
        Line-search sweeps per stage (>= 1).
    parallel : bool
        Run stages 2+ with a thread pool of `max_workers`. Stage 1 is chained
        and always sequential.
    max_workers : int
        Number of threads when `parallel` is True.
    backend : str
        "fsl_mcflirt" (FSL optimiser port, fastest), "sitk_mcflirt" (ITK Powell
        mimicking it).
    verbose : bool
        Print the scaled schedule and warnings.
    only_print_schedule: bool
        Utility to simply print the estimated schedule parameters without conducting
        any registration.

    Returns
    -------
    list of sitk.Euler3DTransform
        One per volume, in ITK convention: each maps reference-space points to
        the corresponding volume's space.
    """
    if not isinstance(sweeps, int) or sweeps < 1:
        raise ValueError(f"sweeps must be an integer >= 1, got {sweeps}")

    if brain_scale is None:
        if brain_size_mask:
            brain_scale = brain_size_mm(brain_size_mask) / HUMAN_BRAIN_SIZE_REF
            if verbose:
                print(f"Parameters are scaled by {brain_scale:.4f}X relative to human brain size.")
        else:
            brain_scale = 1
            if verbose:
                print("Warning: no brain_scale or brain_size_mask given; using human-scale "
                      "parameters (brain_scale = 1).")

    schedule = _read_schedule(schedule)

    # scale lengths (spacing, translations) but not angles (rotations);
    # tol is ordered like the Euler3D parameters: rx, ry, rz, tx, ty, tz
    schedule_scaled = [
        (spacing * brain_scale, np.array([rot_tol] * 3 + [trans_tol * brain_scale] * 3))
        for (spacing, rot_tol, trans_tol) in schedule
    ]

    if verbose:
        print(f"Schedule (brain_scale = {brain_scale:.4f})")
        print(f"{'stage':>5} | {'spacing (mm)':>12} | {'rot tol (deg)':>13} | {'trans tol (mm)':>14}")
        print("-" * 54)
        for stage, (spacing, tol) in enumerate(schedule_scaled):
            print(f"{stage + 1:>5} | {spacing:>12.4f} | {np.rad2deg(tol[0]):>13.4f} | {tol[3]:>14.5f}")

    if only_print_schedule:
        return schedule_scaled

    if backend == "sitk_mcflirt":
        _register_pair = functools.partial(register_pair_sitk_mcflirt, metric=metric, sweeps=sweeps)
    elif backend == "fsl_mcflirt":
        _register_pair = functools.partial(register_pair_fsl_mcflirt, metric=metric, sweeps=sweeps)
    else:
        raise ValueError(f"unknown backend {backend}")

    moving_img = _as_float_image(moving_img)
    ref_img = _as_float_image(ref_img)
    if ref_img.GetDimension() != 3:
        raise ValueError(f"ref_img input must be 3D, got {ref_img.GetDimension()}D instead")
    if moving_img.GetDimension() != 4:
        raise ValueError(f"moving_img input must be 4D, got {moving_img.GetDimension()}D instead")

    size_4d = moving_img.GetSize()
    num_volumes = size_4d[3]
    volumes = []
    for i in range(num_volumes):
        extractor = sitk.ExtractImageFilter()
        extractor.SetSize([size_4d[0], size_4d[1], size_4d[2], 0])
        extractor.SetIndex([0, 0, 0, i])
        volumes.append(extractor.Execute(moving_img))

    center = sitk.CenteredTransformInitializer(
        ref_img, ref_img, sitk.Euler3DTransform(),
        sitk.CenteredTransformInitializerFilter.MOMENTS,
    ).GetCenter()
    identity = sitk.Euler3DTransform()
    identity.SetCenter(center)

    mid = num_volumes // 2
    transforms = [None] * num_volumes

    for stage, (spacing, tol) in enumerate(schedule_scaled):
        fixed = isotropic_upsample_and_pad(
            ref_img, interpolation=sitk.sitkBSpline5, clip_negative=True, spacing=spacing, pad=padding
        )

        # slice count of the stage grid itself, excluding padded edge slices
        if verbose and fixed.GetSize()[2] - 2 * padding < 3:
            print(f"Warning: stage {stage + 1} grid has < 3 slices "
                  "(MCFLIRT would switch to 2D mode here)")

        if stage == 0:
            # chained: middle volume from identity, then outward
            transforms[mid] = _register_pair(fixed, volumes[mid], identity, tol)
            with tqdm(total=num_volumes - 1, desc=f"stage {stage + 1}") as pbar:
                for order in (range(mid + 1, num_volumes), range(mid - 1, -1, -1)):
                    prev = transforms[mid]
                    for i in order:
                        transforms[i] = _register_pair(fixed, volumes[i], prev, tol)
                        prev = transforms[i]
                        pbar.update(1)
        elif parallel:
            # each volume from its own previous estimate -> parallel
            with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
                with tqdm(total=num_volumes, desc=f"stage {stage + 1}") as pbar:
                    futures = {
                        executor.submit(_register_pair, fixed, volumes[i], transforms[i], tol): i
                        for i in range(num_volumes)
                    }
                    for future in concurrent.futures.as_completed(futures):
                        transforms[futures[future]] = future.result()
                        pbar.update(1)
        else:
            for i in tqdm(range(num_volumes), desc=f"stage {stage + 1}"):
                transforms[i] = _register_pair(fixed, volumes[i], transforms[i], tol)

    return transforms


# ----------------------------------------------------------------------------
# Alternative backend: port of FSL MISCMATHS optimise.cc (used by MCFLIRT)
# ----------------------------------------------------------------------------
# All code below this banner was written by Claude (Anthropic's AI assistant),
# September 2026, as a Python port of FSL source code:
#   * _estquadmin, _extrapolatept, _nextpt, _findinitialbound and _optimise1d
#     are line-by-line ports of the functions of the same names in
#     miscmaths/optimise.cc, from the FSL 5.0.1 source distribution
#     (FMRIB, University of Oxford), retrieved from the FreeSurfer-hosted
#     mirror: https://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/
#     tutorial_packages_centos6/OSX/fsl_501/src/miscmaths/optimise.cc
#     The only addition is the step guard in _findinitialbound.
#   * How these are driven (one sweep per stage, boundguess = 10 then 1,
#     tolerances, 8 / 4 / 4 mm stages, stage-1 chaining) follows mcflirt.cc
#     from McFLIRT v2.0 (Peter Bannister, FMRIB) and mcflirt/Globaloptions.h
#     from the same FSL 5.0.1 source tree.
#   * _normcorr_cost is not a port: FSL's normcorr implementation
#     (newimage/costfns.cc) was not consulted. It is a numpy 1 - Pearson r
#     over overlapping voxels, with the same optimum.
# FSL is distributed under the FSL licence (free for non-commercial use);
# check its terms before redistributing this derived code.
# ----------------------------------------------------------------------------
def _estquadmin(x1, xm, x2, y1, ym, y2):
    ad = (xm - x2) * (ym - y1) - (xm - x1) * (ym - y2)
    bd = -(xm * xm - x2 * x2) * (ym - y1) + (xm * xm - x1 * x1) * (ym - y2)
    det = (xm - x2) * (x2 - x1) * (x1 - xm)
    if abs(det) > 1e-15 and ad / det < 0:  # quadratic only has a maximum
        return False, 0.0
    if abs(ad) > 1e-15:
        return True, -bd / (2 * ad)
    return False, 0.0


def _extrapolatept(x1, xm, x2):
    r = 0.3819660
    if abs(x2 - xm) > abs(x1 - xm):
        return r * x2 + (1 - r) * xm
    return r * x1 + (1 - r) * xm


def _nextpt(x1, xm, x2, y1, ym, y2):
    ok, xn = _estquadmin(x1, xm, x2, y1, ym, y2)
    if (not ok) or xn < min(x1, x2) or xn > max(x1, x2):
        xn = _extrapolatept(x1, xm, x2)
    return xn


def _findinitialbound(x1, xm, y1, ym, f, max_steps=100):
    ef, maxextrap = 1.6, 3.2
    if y1 < ym:
        x1, xm, y1, ym = xm, x1, ym, y1
    d = -1.0 if xm < x1 else 1.0
    x2 = xm + ef * (xm - x1)
    y2 = f(x2)
    steps = 0
    while ym > y2 and steps < max_steps:  # step guard is ours, not FSL's
        steps += 1
        maxx2 = xm + maxextrap * (x2 - xm)
        ok, nx = _estquadmin(x1, xm, x2, y1, ym, y2)
        if (not ok) or (nx - x1) * d < 0 or (nx - maxx2) * d > 0:
            nx = xm + ef * (x2 - x1)
        ny = f(nx)
        if (nx - xm) * (nx - x1) < 0:  # nx between x1 and xm
            if ny < ym:
                x2, y2, xm, ym = xm, ym, nx, ny
                break
            x1, y1 = nx, ny
        else:  # nx between xm and maxx2
            if ny > ym:
                x2, y2 = nx, ny
                break
            elif (nx - x2) * d < 0:
                x1, y1, xm, ym = xm, ym, nx, ny
            else:
                x1, y1, xm, ym, x2, y2 = xm, ym, x2, y2, nx, ny
    return x1, xm, x2, y1, ym, y2


def _optimise1d(f, unittol, y0, boundguess=10.0, max_iter=100):
    xm, ym = 0.0, y0
    x1 = boundguess * unittol
    y1 = f(x1)
    x1, xm, x2, y1, ym, y2 = _findinitialbound(x1, xm, y1, ym, f)
    min_dist = 0.1 * unittol
    it = 0
    while True:
        it += 1
        if it > max_iter or abs((x2 - x1) / unittol) <= 1.0:
            break
        xn = _nextpt(x1, xm, x2, y1, ym, y2)
        dn = -1.0 if x2 < x1 else 1.0
        if abs(xn - x1) < min_dist:
            xn = x1 + dn * min_dist
        if abs(xn - x2) < min_dist:
            xn = x2 - dn * min_dist
        if abs(xn - xm) < min_dist:
            xn = _extrapolatept(x1, xm, x2)
        if abs(xm - x1) < 0.4 * unittol:
            xn = xm + dn * 0.5 * unittol
        if abs(xm - x2) < 0.4 * unittol:
            xn = xm - dn * 0.5 * unittol
        yn = f(xn)
        if (xn - xm) * (x2 - xm) > 0:  # keep xn between x1 and xm
            x1, x2, y1, y2 = x2, x1, y2, y1
        if yn < ym:
            x2, y2, xm, ym = xm, ym, xn, yn
        else:
            x1, y1 = xn, yn
    return xm, ym


def _normcorr_cost(fixed, moving, center):
    """1 - Pearson r over fixed voxels that map inside the moving FOV."""
    ref = sitk.GetArrayFromImage(fixed).ravel().astype(np.float64)
    tf = sitk.Euler3DTransform()
    tf.SetCenter(center)

    def cost(p):
        tf.SetParameters(tuple(float(v) for v in p))
        warped = sitk.Resample(moving, fixed, tf, sitk.sitkLinear,
                               float("nan"), sitk.sitkFloat32)
        m = sitk.GetArrayFromImage(warped).ravel()
        ok = np.isfinite(m)
        if ok.sum() < 10:
            return 1.0
        r = np.corrcoef(ref[ok], m[ok])[0, 1]
        return 1.0 - r if np.isfinite(r) else 1.0

    return cost