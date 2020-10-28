# check https://fred.stlouisfed.org/docs/api/fred/series_observations.html for troubleshooting API calls

if (!require("rjson")) install.packages("rjson")
if (!require("dplyr")) install.packages("dplyr")

r_fred <- 
  function( series , 
            fred_key = Sys.getenv("FRED_API_KEY") , # put FRED API key here, defaults to looking for an environment variable
            start_date = "1776-07-04", # default
            last_date = "9999-12-31", # default
            freq = "m" , # see link at top for list of frequencies
            agg_mthd = "avg"){
    require(rjson)
    require(dplyr)
    
    query = paste0("https://api.stlouisfed.org/fred/series/observations?",
                   paste("series_id", sep="=", series) ,
                   sep = "&" ,
                   paste("observation_start" , sep = "=", start_date ) ,
                   sep = "&" ,
                   paste("observation_end", sep = "=" , last_date) ,
                   sep = "&" ,
                   paste("frequency", sep = "=" , freq) ,
                   sep = "&" ,
                   paste("aggregation_method", sep = "=" , agg_mthd) ,
                   sep = "&" ,
                   paste("api_key" ,sep="=", fred_key) ,
                   sep = "&" ,
                   "file_type=json"
    )
    
    response_json = httr::GET(url = query) %>% httr::content("text")
    
    response = fromJSON( response_json )
    
    output = bind_rows(response$observations) %>% select( date, value) %>% mutate(value=as.numeric(value))
    
    return(output)
    
  }
