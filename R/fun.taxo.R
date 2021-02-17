#!/usr/bin/Rscript

# ==============================================================================
# authors         :Ghislain Vieilledent, Aurelien Colas
# email           :ghislain.vieilledent@cirad.fr, aurelien.colas@insa-lyon.fr
# license         :GPLv3
# ==============================================================================

# ==================
# Taxonomy
# ==================

fun.taxo <- function(path,name,spdir,spname,enough,npix){

  ## Taxonomy
  tax.data <- data.frame(binomial=NA,authority=NA,kingdom=NA,family=NA,iucn=NA)
  tax.data$binomial <- spname
  eolid.df <- tryCatch(get_eolid_(spname)[[1]],error = function(e){NA})
  if(is.null(eolid.df)||is.na(eolid.df)){
    tax.data$kingdom <- NA
    tax.data$authority <- NA
    tax.data$family <- NA
    tax.data$iucn <- NA
    text.cut <- "The species was not found in the Encyclopedia of Life (EOL). More information on EOL's website at http://eol.org."
    jpeg(file=paste0(path,"/imagesquare.jpg"),width=100,height=100,units="px")
    grid.raster(matrix(rep(grey(0.9),100), ncol=10),interpolate=FALSE)
    dev.off()
  } else {
    if ("EOL Dynamic Hierarchy" %in% eolid.df$source) {
      eol.id <- eolid.df %>% filter(source=="EOL Dynamic Hierarchy") %>% pull(eolid) %>% first()
      eolpage.id <- eolid.df %>% filter(source=="EOL Dynamic Hierarchy") %>% pull(pageid) %>% first()
    } else if ("IUCN" %in% eolid.df$source) {
      eol.id <- eolid.df %>% filter(source=="IUCN") %>% pull(eolid) %>% first()
      eolpage.id <- eolid.df %>% filter(source=="IUCN") %>% pull(pageid) %>% first()
    } else if ("NCBI" %in% eolid.df$source) {
      eol.id <- eolid.df %>% filter(source=="NCBI") %>% pull(eolid) %>% first()
      eolpage.id <- eolid.df %>% filter(source=="NCBI") %>% pull(pageid) %>% first()
    } else  {
      eol.id <- eolid.df %>% pull(eolid) %>% first()
      eolpage.id <- eolid.df %>% pull(pageid) %>% first()
    }
    eolclassif.df <- classification(eol.id,db="eol")[[1]]

    if(!(is.character(eolclassif.df))){
      if(nrow(eolclassif.df[(eolclassif.df$rank=="kingdom")&!(is.na(eolclassif.df$rank)),])>0){tax.data$kingdom <- eolclassif.df %>% filter(rank=="kingdom") %>% pull(name)}
      if(nrow(eolclassif.df[(eolclassif.df$rank=="family")&!(is.na(eolclassif.df$rank)),])>0){tax.data$family <- eolclassif.df %>% filter(rank=="family") %>% pull(name)}
      ## Authority
      if (tax.data$kingdom %in% c("Animalia","Metazoa")) {
        tax.data$kingdom <- "Animalia"
        gnr <- gnr_resolve(spname,data_source_ids=3)
        if(nrow(gnr)>0){
          gnr <- gnr %>% pull(matched_name)
          matched_name <- regmatches(gnr,regexpr("[[:alpha:]]+ [[:alpha:]]+", gnr))[1]
          auth <- regmatches(gnr,regexpr("\\([[:alpha:]]+[,[:space:]]+[[:digit:]]{0,4}\\)", gnr))
          authority <- ifelse(length(auth)==0,NA,auth)
        } else { authority <- NA }
      } else if (tax.data$kingdom %in% c("Plantae","Viridiplantae")) {
        tax.data$kingdom <- "Plantae"
        tnrs_query <- tnrs(query=spname,source="iPlant_TNRS")
        matched_name <- tnrs_query$matchedname
        auth <- tnrs_query$authority
        authority <- ifelse(length(auth)==0,NA,auth)
      } else(authority <- NA )
      if(!(is.na(authority))){
        tax.data$authority <- as.character(ifelse(matched_name==spname,iconv(iconv(authority,from="UTF-8",to="ASCII//TRANSLIT"),from="ASCII//TRANSLIT",to="UTF-8"),NA))
      }
    }
    ## IUCN conservation status
    iucn.summary <- iucn_summary(spname)[[1]]
    tax.data$iucn <- ifelse(!is.na(iucn.summary[1]),as.character(iucn.summary$status),NA)

    #Text alternative
    textFile <- paste0("../figures/",name,"/",spdir,"/desc.txt")
    if(file.exists(textFile)){
      text <- readChar(textFile, file.info(textFile)$size)
      id.pts <- unlist(gregexpr(pattern="\\.(\\s[[:upper:]]|$)",text))
      id.cut <- ifelse(id.pts[1]==-1,1000,max(id.pts[id.pts<=1000]))
      text.cut <- substring(text,1,id.cut)
      text.cut <- paste0(text.cut," $[\\dots]$"," https://lemursofmadagascar.com")
    } else {
      text.cut <- paste0("Description text was not found in the Encyclopedia of Life (EOL). More information on EOL's website at http://eol.org/",eolpage.id)
    }
    # Text
    # text.id <- eol_pages(eolpage.id,texts_page=1,vetted=2, detail=T)$data_objects
    # if(is.data.frame(text.id)){
    #   curl::curl_download(url=paste0("eol.org/data_objects/",text.id$dataobjectversionid),destfile=paste0(path,"/text.html"))
    #   # Convert HTML to text
    #   html <- htm2txt(iconv(iconv(readLines(paste0(path,"/text.html")),from="UTF-8",to="ASCII//TRANSLIT"),from="ASCII//TRANSLIT",to="UTF-8"))
    #   i=1
    #   while (nchar(html[i])<90) {
    #     i=i+1
    #   }
    #   if((i<66)&(substr(html[i],1,1)!="?")){
    #     text <- html[i]
    #     # We select the complete sentences with less than 1000 symbols. We also add dots and link to eol page
    #     id.pts <- unlist(gregexpr(pattern="\\.(\\s[[:upper:]]|$)",text))
    #     id.cut <- ifelse(id.pts[1]==-1,1000,max(id.pts[id.pts<=1000]))
    #     text.cut <- substring(text,1,id.cut)
    #     if(!(is.na(text.cut))){
    #       text.cut <- paste0(text.cut," $[\\dots]$"," http://eol.org/",eolpage.id)
    #     } else {
    #       text.cut <- paste0("Description text was not found in the Encyclopedia of Life (EOL). More information on EOL's website at http://eol.org/",eolpage.id)
    #     }
    #   }else {
    #     text.cut <- paste0("Description text was not found in the Encyclopedia of Life (EOL). More information on EOL's website at http://eol.org/",eolpage.id)
    #   }
    # } else {
    #   text.cut <- paste0("Description text was not found in the Encyclopedia of Life (EOL). More information on EOL's website at http://eol.org/",eolpage.id)
    # }

    #Image alternative

    imgFile <- paste0("../figures/",name,"/",spdir,"/",gsub(" ","-",spname),".jpg")
    if(file.exists(imgFile)){
      img <- magick::image_read(imgFile)
      img.bmp <- readbitmap::read.bitmap(imgFile)
      h <- dim(img.bmp)[1]
      w <- dim(img.bmp)[2]
      if (w>=h) {
        img.square <- magick::image_crop(img, paste0(h,"x",h,"+",(w-h)/2))
      } else {
        img.square <- magick::image_crop(img, paste0(w,"x",w,"+","0+",(h-w)/2))
      }
      magick::image_write(img.square,path=paste0(path,"/imagesquare.jpg"),format="jpg")
    } else {
      jpeg(file=paste0(path,"/imagesquare.jpg"),width=100,height=100,units="px")
      grid.raster(matrix(rep(grey(0.9),100), ncol=10),interpolate=FALSE)
      dev.off()
    }

    # # Image
    # img.id <- eol_pages(eolpage.id,images_page=10,vetted=2)$data_objects
    # if(is.data.frame(img.id)){
    #   curl::curl_download(url=paste0("eol.org/data_objects/",img.id$dataobjectversionid),destfile=paste0(path,"/img.html"))
    #   # Taking the lines of the HTML file
    #   html <- readLines(paste0(path,"/img.html"))
    #   i=1
    #   while (substring(html[i],1,44)!="<meta content='http://media.eol.org/content/") {
    #     i=i+1
    #   }
    #   # We select the url of the image
    #   id.ext <- unlist(gregexpr(pattern="\\.jpg",html[i]))
    #   img.url <- substring(html[i],16,id.ext+3)
    #   HTTP <- tryCatch(curl::curl_download(url=img.url,destfile=paste0(path,"/image.jpg")),error = function(e){NA})
    #   if(is.na(HTTP)){
    #     jpeg(file=paste0(path,"/imagesquare.jpg"),width=100,height=100,units="px")
    #     grid.raster(matrix(rep(grey(0.9),100), ncol=10),interpolate=FALSE)
    #     dev.off()
    #   } else {
    #     img.path <- paste0(path,"/image.jpg")
    #     img <- magick::image_read(img.path)
    #     img.bmp <- readbitmap::read.bitmap(img.path)
    #     h <- dim(img.bmp)[1]
    #     w <- dim(img.bmp)[2]
    #     if (w>=h) {
    #       img.square <- magick::image_crop(img, paste0(h,"x",h,"+",(w-h)/2))
    #     } else {
    #       img.square <- magick::image_crop(img, paste0(w,"x",w,"+","0+",(h-w)/2))
    #     }
    #     magick::image_write(img.square,path=paste0(path,"/imagesquare.jpg"),format="jpg")
    #   }
    # } else {
    #   jpeg(file=paste0(path,"/imagesquare.jpg"),width=100,height=100,units="px")
    #   grid.raster(matrix(rep(grey(0.9),100), ncol=10),interpolate=FALSE)
    #   dev.off()
    # }
  }

  ##========================================
  ## Save objects to be loaded by knitr
  if(enough){
    save(list=c("tax.data","text.cut"),file=paste0(path,"/taxonomy.rda"))
  } else {
    save(list=c("tax.data","text.cut","npix"),file=paste0(path,"/data.rda"))
  }
}
