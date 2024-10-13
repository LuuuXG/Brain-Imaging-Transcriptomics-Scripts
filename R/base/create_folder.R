#require('crayon')

# Function to create a folder
create_folder <- function(folder_path) {
  if (!dir.exists(folder_path)) {
    dir.create(folder_path)
    cat(note(paste("Folder", folder_path, "has been created.\n")))
  } else {
    cat(warn(paste("Folder", folder_path, "already exists.\n")))
  }
}