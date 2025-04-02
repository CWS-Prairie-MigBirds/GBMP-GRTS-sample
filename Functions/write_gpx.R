# custom function to write GPX file
#adapted by Josiah from https://rdrr.io/github/s-u/snippets/src/R/gpx.R 

write_gpx <- function(lat, lon, time = NULL, name = NULL, out_file) {
  o <- c('<gpx version="1.1" creator="R">')
  
  if (is.null(time)) {
    for (i in seq_along(lat)) {
      if (is.null(name)) {
        o <- c(o, paste('<wpt lat="', lat[i], '" lon="', lon[i], '" />', sep=''))
      } else {
        o <- c(o, paste('<wpt lat="', lat[i], '" lon="', lon[i], '"><name>', name[i], '</name></wpt>', sep=''))
      }
    }
  } else {
    for (i in seq_along(lat)) {
      if (is.null(name)) {
        o <- c(o, paste('<wpt lat="', lat[i], '" lon="', lon[i], '"><time>', 
                        paste(gsub(' ','T', as.character(time[i])), 'Z', sep=''), 
                        '</time></wpt>', sep=''))
      } else {
        o <- c(o, paste('<wpt lat="', lat[i], '" lon="', lon[i], '"><time>', 
                        paste(gsub(' ','T', as.character(time[i])), 'Z', sep=''), '</time><name>', 
                        name[i], '</name></wpt>', sep=''))
      }
    }
  }
  
  o <- c(o, '</gpx>')
  
  if (is.character(out_file) || inherits(out_file, "connection")) 
    cat(o, file=out_file, sep='\n')
}