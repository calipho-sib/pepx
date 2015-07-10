pepx
====

Peptide matcher across nextprot entries and variants

## To build the application

gcc pep.c -o pepx

either generate the indexes or copy them: 'cp /share/sib/common/Calipho/alain/pepxvaridx/* .'

## When running the application

Pepx possible arguments are:

- --search (short=-s) to perform a peptide search
- --build (short=-b) to build indexes (requires the isoform-file name as mandatory argument)
- --version (short=-v) to show current version
- --help (short=-h) to show this help
- --index-folder (short=-x) to specify an index folder (default is .)
- --variant-folder (short=-w) to specify a folder for json variants (required for build command when ignore-variants flag is not set)
- --rest-url (short=-r) REST server to retrieve json variants when variant folder is empty (for 1rst build with a given variant folder)
- --peptide-file (short=-p) a file with peptides to search (1 peptide/line, if not provided peptides will be read from stdin)
- --ignore-variants to build indexes not considering variants
- --IL to build indexes merging I and L
- --noiso (short=-n) to output search results at the entry level

Current limitation:

- poly-AA stretches > 6 cannot be found
- only snp-style (1 AA for 1 other AA), and 1-AA-miss variants are accounted
- only 32 variant accounted within a given x-mer
- only 1 joker (X) allowed in a given x-mer

