I will configure the repository for deployment to GitHub Pages.

# 1. Configuration
- Create a `.gitignore` file to exclude large data files (genomes, archives) and temporary files, ensuring only the tool code and block data are committed.
- Create a GitHub Actions workflow `.github/workflows/deploy.yml` that automatically deploys the contents of the `web/` folder to the `gh-pages` branch whenever code is pushed.

# 2. File Preparation
- Add a `.gitkeep` file to `web/data/pombase_dy47073/` so the empty directory is tracked by git.

# 3. Git Operations
- Create and switch to a new branch `feature/web-deployment`.
- Stage the `web/` directory, `cli/` scripts, and configuration files.
- Commit the changes.
- Push the branch to the remote repository.

# 4. Post-Deployment Instructions
After I complete these steps, the code will be on GitHub. You will just need to:
1.  Go to your repository settings on GitHub.
2.  Navigate to **Pages**.
3.  Ensure the source is set to the `gh-pages` branch (this branch will be created automatically by the action).
