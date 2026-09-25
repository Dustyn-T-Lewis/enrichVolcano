# Clean a pathway name for display

Hallmark and GO slim terms take their label from a table written by hand
and shipped with the package (`inst/extdata/term_labels.csv`): `"short"`
fits two ring lines, `"clean"` spells every word out. The package's
earlier cleaning rules informed the short names. Any other name is
cleaned plainly: the database prefix and any `>` or `__` levels above
the last are dropped, underscores become spaces, words are title-cased
and common acronyms such as DNA and mRNA restored.
[`list_labels()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/list_labels.md)
lists the labels a figure will draw, to edit and pass back through
`labels`.

## Usage

``` r
clean_label(name, width = 15, style = c("short", "clean"))
```

## Arguments

- name:

  Character vector of term names.

- width:

  Wrap width in characters. A label that already holds a line break is
  not wrapped again.

- style:

  `"short"` or `"clean"`. Only table terms differ between the two.

## Value

Character vector of labels.

## Examples

``` r
clean_label(c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "REACTOME_TCA_CYCLE"))
#> [1] "OXPHOS"    "TCA Cycle"
clean_label("HALLMARK_OXIDATIVE_PHOSPHORYLATION", style = "clean", width = 40)
#> [1] "Oxidative Phosphorylation"
```
