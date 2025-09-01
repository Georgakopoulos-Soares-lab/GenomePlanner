# GenomePlanner

**GenomePlanner** is a high-performance C++ tool that uses **dynamic programming** to find the most cost-effective plan for constructing a target DNA sequence using fragments from a source sequence.  

It models the trade-offs between two primary DNA acquisition methods:

- **Replication (PCR):** Reusing an existing DNA fragment found in the source sequence for a fixed cost.  
- **Synthesis:** Building a DNA fragment from scratch at a specified cost per base.  

The tool leverages a **powerful FM-index** (powered by the [SDSL-lite](https://github.com/simongog/sdsl-lite) library) for instantaneous lookups of source fragments, enabling an efficient DP-based search for the globally optimal (i.e., cheapest) construction plan.  

---

## Core Concepts

- **FM-Index:**  
  A compressed full-text index that allows for extremely fast counting and locating of any substring within the source genome. This answers the question: *"Does this piece of DNA exist in my source?"*

- **Dynamic Programming:**  
  The algorithm builds the optimal plan incrementally. It calculates the cheapest possible cost to construct the target sequence up to every single base, ensuring the final plan is the most cost-effective solution possible.

- **Cost Model:**  
  The planner's decisions are entirely driven by a user-defined cost model that balances:
  - The fixed cost of replication (`--pcr`)
  - The per-base cost of synthesis (`--synth`)
  - The cost of joining two DNA fragments together (`--join`)

---

## Prerequisites

To compile and run the planner, you will need:

- A **C++17 compliant compiler** (e.g., `g++`)
- The [SDSL-lite](https://github.com/simongog/sdsl-lite) library
- The **zlib** compression library (`-lz`)

---

## Compilation

Assuming dependencies (like SDSL-lite) are installed in your environment (e.g., via Conda):

```bash
g++ -O3 -std=gnu++17 dp_cost.cpp -o dp_cost     -lsdsl -ldivsufsort -ldivsufsort64 -lz
```

---

## Usage

The program is controlled via command-line arguments:

```bash
./dp_cost --source <SRC.fasta.gz>           --target <TGT.fasta.gz>           --W <INT>           --pcr <FLOAT>           --join <FLOAT>           --synth <FLOAT>           --summary-out <path.tsv>
```

### Flags

| Flag            | Description                                                                 |
|-----------------|-----------------------------------------------------------------------------|
| `--source`      | Path to the source genome FASTA (plain or `.gz`).                           |
| `--target`      | Path to the target genome FASTA (plain or `.gz`).                           |
| `--W`           | The maximum block size (k-mer length) to consider.                          |
| `--pcr`         | The fixed cost for a replication (reuse) operation.                         |
| `--join`        | The cost to join two consecutive blocks.                                    |
| `--synth`       | The per-base cost for a synthesis operation.                                |
| `--summary-out` | **Required.** The path for the output summary TSV file.                     |
| `--moves-out`   | Optional. Path for a detailed step-by-step plan.                            |
| `--progress`    | Optional. Percentage step for progress updates (default: `5`).              |

---

## Example Walkthrough: Favoring Short K-mer Replication

This example demonstrates how to configure the cost model to favor the reuse of very short k-mers (e.g., 8–9 base pairs), even if it results in a highly fragmented plan.

**Source:** `NC_010323.1.fasta.gz` (Rhododendron virus A)  
**Target:** `NC_009899.1.fasta.gz` (Amasya cherry disease associated virus)  

### Run Command

The costs are set so that reusing an 8-mer (`pcr` cost of `3.0` + `join` cost of `0.1` = `3.1`) is cheaper than synthesizing it (`8 bases * synth cost of 0.5 = 4.0`).  
The minimal join cost removes the penalty for fragmentation.

```bash
./dp_cost --source NC_010323.1.fasta.gz           --target NC_009899.1.fasta.gz           --W 200           --pcr 3           --join 0.1           --synth 0.5           --summary-out summary_low_join.tsv
```

---

## Results (`summary_low_join.tsv`)

| Metric                  | Value   |
|-------------------------|---------|
| Chromosome              | NC_009899.1_Paramecium_bursaria_Chlorella_virus_AR158__complete_genome |
| Total Length            | 344,686 |
| Total Cost              | 122,795 |
| Total Blocks            | 38,070  |
| Replication Blocks      | 35,326  |
| Synthesis Blocks        | 2,744   |
| Avg. Replication Length | 9.02 bp |
| Avg. Synthesis Length   | 9.48 bp |

---

## Analysis

The results show a resounding success. The plan consisted of **35,326 replicated blocks** with an **average length of just 9.0 base pairs**.  

This proves that the algorithm correctly identified and utilized the short, random k-mer matches that exist between the two unrelated viruses—just as intended by the cost model.
