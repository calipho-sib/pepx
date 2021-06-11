# neXtProt - The knowledge resource on human proteins

This  is a code repository for the SIB - Swiss Institute of Bioinformatics CALIPHO group neXtProt project

See: http://www.nextprot.org/

# pepx

Peptide matcher across neXtProt entries and variants

## To build the application

gcc pepx.c -o pepx

indexing exemple:  pepx --IL --ignore-variants -x my_index_folder -b my_sequence_file
searching exemple: pepx --IL -x my_index_folder -s NFIVSTWHR

## When running the application

Pepx possible arguments are:

- --search (short=-s) to perform a peptide search
- --build (short=-b) to build indexes (requires the isoform-file name as mandatory argument)
- --version (short=-v) to show current version
- --help (short=-h) to show this help
- --index-folder (short=-x) to specify an index folder (default is .)
- --variant-folder (short=-w) to specify a folder for json variants (required for build command when ignore-variants flag is not set)
- --peptide-file (short=-p) a file with peptides to search (1 peptide/line, if not provided peptides will be read from stdin)
- --ignore-variants to build indexes not considering variants
- --IL to build indexes merging I and L and search these indexes
- --7mers-only to build only 7-mer indexes (saves disk space)
- --noiso (short=-n) to output search results at the entry level

Current limitation:

- poly-AA stretches > 7 can not be exactly found and are overpredicted
- only snp-style (1 AA for 1 other AA), and 1-AA-miss variants are accounted
- max 128 variants accounted within a given x-mer
- only 1 joker (X) allowed in a given x-mer in search mode

