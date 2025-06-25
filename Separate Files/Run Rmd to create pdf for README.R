# Run the Rmd and md file

rmarkdown::render(
  "C:/users/u1107832/OneDrive - Australian National University/R Package Development/groupedHG/Separate Files/README_pdf_document.Rmd",
  output_file = "README_pdf_document.pdf",
  output_dir = "C:/users/u1107832/OneDrive - Australian National University/R Package Development/groupedHG"
)


rmarkdown::render(
  "C:/users/u1107832/OneDrive - Australian National University/R Package Development/groupedHG/Separate Files/HG-Vignette_pdf_document.Rmd",
  output_file = "HG_Vignette_pdf_document.pdf",
  output_dir = "C:/users/u1107832/OneDrive - Australian National University/R Package Development/groupedHG"
)
