require(httr)
require(jsonlite)

url <- "https://biotools.fr/human/ensembl_symbol_converter/"
ids <- scan("HPA-genes-ENSEMBL", character(), quote = "")
ids_json <- toJSON(ids)

body <- list(api=1, ids=ids_json)
r <- POST(url, body = body)

output <- fromJSON( content(r, "text"), flatten=TRUE)
outfile <- file("HPA-gene-symbols")
writeLines(unlist(output, use.names=FALSE), outfile)
close(outfile)
