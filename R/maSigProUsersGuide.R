maSigProUsersGuide<-function (view = TRUE) 
{
    f <- system.file("doc", "maSigProUsersGuide.pdf", package = "maSigPro")
    if (view) {
        if (.Platform$OS.type == "windows") 
            shell.exec(f)
        else system(paste(Sys.getenv("R_PDFVIEWER"), f, "&"))
    }
    return(f)
}

