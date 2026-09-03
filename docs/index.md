```{include} ../README.md
```

## Finding your way around this documentation

This documentation is organised around what you are trying to do right now.

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item-card} {octicon}`rocket;1.5em;sd-mr-1` Tutorial
:link: tutorials/index
:link-type: doc

**Start here if RABIES is new to you.** A guided run of the complete
pipeline on a small example dataset, from raw BIDS input to a
connectivity map.
:::

:::{grid-item-card} {octicon}`tools;1.5em;sd-mr-1` How-to guides
:link: how_to/index
:link-type: doc

**Practical recipes for a specific goal.** Installing RABIES, running it
in a container, tuning a failed registration, designing a confound
correction strategy, contributing code.
:::

:::{grid-item-card} {octicon}`book;1.5em;sd-mr-1` Reference
:link: reference/index
:link-type: doc

**Look something up.** Command line options, input requirements, every
output file, and the precise definition of every metric.
:::

:::{grid-item-card} {octicon}`light-bulb;1.5em;sd-mr-1` Explanation
:link: explanation/index
:link-type: doc

**Understand how RABIES works and why.** The preprocessing and confound
correction workflows, the analyses, and the data quality framework.
:::

::::

```{toctree}
---
maxdepth: 2
hidden:
caption: Getting started
---
tutorials/index
```

```{toctree}
---
maxdepth: 2
hidden:
caption: How-to guides
---
how_to/index
```

```{toctree}
---
maxdepth: 2
hidden:
caption: Reference
---
reference/index
```

```{toctree}
---
maxdepth: 2
hidden:
caption: Explanation
---
explanation/index
```
