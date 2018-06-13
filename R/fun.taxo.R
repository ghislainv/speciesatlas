#!/usr/bin/Rscript

# ==============================================================================
# authors         :Ghislain Vieilledent, Aurélien Colas
# email           :ghislain.vieilledent@cirad.fr, aurelien.colas@insa-lyon.fr
# license         :GPLv3
# ==============================================================================

# ==================
# Taxonomy
# ==================

fun.taxo <- function(path,name,spdir,spname,enough,npix){
  
  ## Taxonomy
  eolid.df <- get_eolid_(spname)[[1]]
  tax.data <- data.frame(binomial=NA,authority=NA,kingdom=NA,family=NA,iucn=NA)
  tax.data$binomial <- spname
  if(is.null(eolid.df)){
    tax.data$kingdom <- NA
    tax.data$authority <- NA
    tax.data$family <- NA
    tax.data$iucn <- NA
    text.cut <- "The species was not found in the Encyclopedia of Life (EOL). More information on EOL's website at http://eol.org."
    jpeg(file=paste0(path,"/image_square.jpg"),width=100,height=100,units="px")
    grid.raster(matrix(rep(grey(0.9),100), ncol=10),interpolate=FALSE)
    dev.off()
  } else {
    if ("IUCN" %in% eolid.df$source) {
      eol.id <- eolid.df %>% filter(source=="IUCN") %>% pull(eolid) %>% first()
      eolpage.id <- eolid.df %>% filter(source=="IUCN") %>% pull(pageid) %>% first()
    } else if ("NCBI" %in% eolid.df$source) {
      eol.id <- eolid.df %>% filter(source=="NCBI") %>% pull(eolid) %>% first()
      eolpage.id <- eolid.df %>% filter(source=="NCBI") %>% pull(pageid) %>% first()
    } else {
      eol.id <- eolid.df %>% pull(eolid) %>% first()
      eolpage.id <- eolid.df %>% pull(pageid) %>% first()
    }
    eolclassif.df <- classification(eol.id,db="eol")[[1]]
    tax.data$kingdom <- eolclassif.df %>% filter(rank=="kingdom") %>% pull(name)
    tax.data$family <- eolclassif.df %>% filter(rank=="family") %>% pull(name)
    ## Authority
    if (tax.data$kingdom %in% c("Animalia","Metazoa")) {
      tax.data$kingdom <- "Animalia"
      gnr <- gnr_resolve(spname,data_source_ids=3) %>% pull(matched_name)
      matched_name <- regmatches(gnr,regexpr("[[:alpha:]]+ [[:alpha:]]+", gnr))
      auth <- regmatches(gnr,regexpr("\\([[:alpha:]]+[,[:space:]]+[[:digit:]]{0,4}\\)", gnr))
      authority <- ifelse(length(auth)==0,NA,auth)
    } else if (tax.data$kingdom %in% c("Plantae","Viridiplantae")) {
      tax.data$kingdom <- "Plantae"
      tnrs_query <- tnrs(query=spname,source="iPlant_TNRS")
      matched_name <- tnrs_query$matchedname
      auth <- tnrs_query$authority
      authority <- ifelse(length(auth)==0,NA,auth)
    }
    tax.data$authority <- as.character(ifelse(matched_name==spname,gsub('ü','u',authority),NA))
    ## IUCN conservation status
    iucn.summary <- iucn_summary(spname)[[1]]
    tax.data$iucn <- ifelse(!is.na(iucn.summary[1]),as.character(iucn.summary$status),NA)
    
    ## Encyclopedia Of Life (EOL)
    eol.page <- Reol::DownloadEOLpages(eolpage.id,to.file=FALSE)
    if (is.null(eol.page)) {
      text.cut <- "The species was not found in the Encyclopedia of Life (EOL). More information on EOL's website at http://eol.org."
      jpeg(file=paste0(path,"/image_square.jpg"),width=100,height=100,units="px")
      grid.raster(matrix(rep(grey(0.9),100), ncol=10),interpolate=FALSE)
      dev.off()
    } else {
      DOI <- Reol::GatherDataObjectInformation(eol.page)
      # Text
      w.text <- grep("Text",DOI$dataType)
      if(!is.na(w.text[1])){
        text <- DOI$description[w.text[1]]
        # Convert HTML to text
        text.txt <- htm2txt::htm2txt(text)
        # We select the complete sentences with less than 1000 symbols. We also add dots and link to eol page
        id.pts <- unlist(gregexpr(pattern="\\.(\\s[[:upper:]]|$)",text.txt))
        id.cut <- ifelse(id.pts[1]==-1,700,max(id.pts[id.pts<=700]))
        text.cut <- paste0(substring(text.txt,1,id.cut)," $[\\dots]$"," http://eol.org/",eolpage.id)
      } else { 
        text.cut <- "Description text was not found in the Encyclopedia of Life (EOL). More information on EOL's website at http://eol.org."
      }
      
      # Image
      w.image <- grep("image/jpeg",DOI$mimeType)
      image.url <- DOI$mediaURL[w.image[1]]
      if(is.null(image.url)){image.url <- NA}
      if(!is.na(image.url)){
        ext <- extension(image.url)
        HTTP <- tryCatch(curl::curl_download(url=image.url,destfile=paste0(path,"/image",ext)),error = function(e){NA})
        if(is.na(HTTP)){
          jpeg(file=paste0(path,"/image_square.jpg"),width=100,height=100,units="px")
          grid.raster(matrix(rep(grey(0.9),100), ncol=10),interpolate=FALSE)
          dev.off()  
        } else {
          img.path <- paste0(path,"/image",ext)
          img <- magick::image_read(img.path)
          img.bmp <- readbitmap::read.bitmap(img.path)
          h <- dim(img.bmp)[1]
          w <- dim(img.bmp)[2] 
          if (w>=h) {
            img.square <- magick::image_crop(img, paste0(h,"x",h,"+",(w-h)/2))
          } else {
            img.square <- magick::image_crop(img, paste0(w,"x",w,"+","0+",(h-w)/2))
          }
          magick::image_write(img.square,path=paste0(path,"/image_square.jpg"),format="jpg")
        }
      } else {
        jpeg(file=paste0(path,"/image_square.jpg"),width=100,height=100,units="px")
        grid.raster(matrix(rep(grey(0.9),100), ncol=10),interpolate=FALSE)
        dev.off()
      }
    }
  }
  
  ##========================================
  ## Save objects to be loaded by knitr
  if(enough){
    save(list=c("tax.data","text.cut"),file=paste0(path,"/taxonomy.rda"))
  } else {
    save(list=c("tax.data","text.cut","npix"),file=paste0(path,"/data.rda"))
  }    
}