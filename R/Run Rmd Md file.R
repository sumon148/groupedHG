# Run the Rmd and md file

rmarkdown::render(
  "C:/users/u1107832/OneDrive - Australian National University/R Package Development/groupedHG/NEWS.md",
  output_file = "NEWS.pdf",
  output_dir = "C:/users/u1107832/OneDrive - Australian National University/R Package Development/groupedHG"
)

rmarkdown::render(
  "C:/users/u1107832/OneDrive - Australian National University/R Package Development/groupedHG/README.md",
  output_file = "README.pdf",
  output_dir = "C:/users/u1107832/OneDrive - Australian National University/R Package Development/groupedHG"
)


rmarkdown::render(
  "C:/users/u1107832/OneDrive - Australian National University/R Package Development/groupedHG/README.Rmd",
  output_file = "README.pdf",
  output_dir = "C:/users/u1107832/OneDrive - Australian National University/R Package Development/groupedHG"
)
