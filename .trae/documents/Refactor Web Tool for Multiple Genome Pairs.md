Refactoring the web tool to support multiple genome pairs dynamically is the most efficient way to host this on GitHub Pages. Instead of duplicating the site for each pair, we will make the single "app" load different data files based on a user selection.

Here is the implementation plan:

# 1. Data Structure Organization
We will create a structured data directory within the `web/` folder to host the block files for different pairs. This keeps the GitHub Pages site self-contained.
- Create `web/data/pombase_leupold/`
- Create `web/data/pombase_dy47073/`
- Move/Copy the existing `A_to_B.blocks.tsv` and `B_to_A.blocks.tsv` (currently in the project root) into `web/data/pombase_leupold/`.
- Create placeholder files or ensure the directory exists for `pombase_dy47073`.

# 2. Update `web/index.html`
- Add a **Genome Pair Selector** (dropdown menu) at the top of the UI.
- Remove hardcoded "PomBase" and "Leupold" text from the headers and buttons where appropriate, or make them update dynamically via JavaScript.

# 3. Refactor `web/app.js`
- **Configuration**: Add a `GENOME_PAIRS` configuration object mapping IDs (e.g., `pombase_leupold`) to their display names and file paths.
- **State Management**: Add a `currentPair` state variable.
- **Dynamic Loading**: Update `loadDefaultAB` and `loadDefaultBA` to accept file paths or read from `currentPair`.
- **UI Updates**: Create a function to update button labels and status text (e.g., "Liftover PomBase→DY47073") when the selection changes.
- **Versioning**: Update `APP_VERSION` to `v1.01`.

# 4. Execution Steps
1.  Create the directory structure `web/data/`.
2.  Copy existing block files to the `pombase_leupold` folder.
3.  Modify `web/index.html` to add the `<select>` element.
4.  Modify `web/app.js` to implement the switching logic and config.
