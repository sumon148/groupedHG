
# Steps for R Package Development

# Step 1: Set Up a Directory
# Create a folder named after your package (e.g., mypackage).
# Inside this folder, place your code, data, and documentation files.

#Step 1:2. Install Development Tools
# Install the devtools, roxygen2, usethis, and testthat packages:

install.packages(c("devtools", "roxygen2", "usethis", "testthat"))

# Step 3. Create the Package Skeleton
# Already done

# Step 4. Write Functions
# Already done
# In R folder only the function will be available
# For running the package use separate file and put them as application folder

# 3 Step 5. Add Metadata
# Edit the DESCRIPTION file:

#  Fill in fields like Package, Title, Version, Author, Maintainer, and Description.
# Specify dependencies using the Imports or Depends fields.

# 6. Document Your Code
# Use roxygen2 to create documentation:
#
#   Add #' comments above each function.
# Already done
# Run devtools::document() to generate .Rd files in the man/ directory.
devtools::document()

# Step 7. Include Tests
# Set up a testing framework with testthat:
usethis::use_testthat()
# Create test files in the tests/testthat/ directory and write test cases:

# Step 8. Add a README
# Generate a README file using usethis:

usethis::use_readme_rmd()

# Step 9. Check Your Package
# Check for issues using:
devtools::check()
# Resolve any warnings or errors.

# Step 10. Build and Install the Package
# Build your package:
  devtools::build()

# Install it locally:
devtools::install()

#





# To start the R package project:

# Restart R and Clear Workspace
.rs.restartR()
rm(list = ls())
gc()

devtools::check()
