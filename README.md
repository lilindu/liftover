# SchizoLift

**SchizoLift** is a client-side web tool for genome coordinate conversion for _Schizosaccharomyces_ species. It allows researchers to easily convert genomic coordinates between different reference genomes and clade consensus sequences.

**Live Tool:** [https://lilindu.github.io/liftover/](https://lilindu.github.io/liftover/)

## Features

*   **Multi-Genome Support:** Easily switch between different genome pairs (e.g., PomBase ↔ Leupold Consensus).
*   **Client-Side Execution:** Runs entirely in your web browser. No data is sent to a server, ensuring privacy and speed.
*   **Detailed Reporting:** Handles simple point conversion and complex interval stitching, reporting gaps and block boundaries.
*   **Round-Trip Validation:** Includes a feature to map coordinates A→B and then B→A to verify consistency.

## Supported Genome Pairs

*   **PomBase ↔ Leupold Consensus**: Mapping between the standard *S. pombe* reference (PomBase) and the Leupold clade consensus.
*   **PomBase ↔ DY47073**: (Coming Soon) Mapping between PomBase and the DY47073 strain.

## Usage

1.  **Select Genome Pair**: Use the dropdown menu at the top to choose the desired genome pair.
2.  **Input Coordinates**: Enter coordinates in the text area. Supported formats (1-based, inclusive):
    *   `CHR:POS` (e.g., `I:1000`)
    *   `CHR:START-END` (e.g., `II:200-500`)
3.  **Run Liftover**: Click the appropriate button (e.g., "Liftover A→B") to perform the conversion.
4.  **View/Download**: Results are displayed in a table below and can be downloaded as a TSV file.

## Development

The project is hosted on GitHub Pages. The core logic resides in the `web/` directory.

### Directory Structure

*   `web/index.html`: The main user interface.
*   `web/app.js`: Application logic for coordinate conversion and UI management.
*   `web/data/`: Contains the block mapping files for different genome pairs.

### Adding New Genome Pairs

To add a new pair (e.g., GenomeX ↔ GenomeY):

1.  Create a new directory in `web/data/` (e.g., `web/data/genomex_genomey/`).
2.  Generate the block files `A_to_B.blocks.tsv` and `B_to_A.blocks.tsv` using the provided CLI tools (see `cli/` folder).
3.  Place the block files in the new directory.
4.  Update `GENOME_PAIRS` in `web/app.js` to register the new pair.

## License

[MIT License](LICENSE)
