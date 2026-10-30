
# Travel Choice Modeling

This repository documents a Python-based project focused on travel behaviour analysis and discrete choice modeling using the **London Passenger Mode Choice (LPMC)** dataset.

The project explores how travel times, costs, and individual and household characteristics relate to travel mode choices. It also provides a framework for estimating and comparing discrete choice models and evaluating their behavioural implications.

## Objectives

- Explore individual travel behaviour and mode choice patterns.
- Estimate discrete choice models using observed travel decisions.
- Examine the influence of travel time, cost, and socioeconomic characteristics on mode choice.
- Compare model performance, statistical properties, and behavioural responses.
- Investigate the potential applications of choice models in transportation planning.

## Dataset

The project uses the **London Passenger Mode Choice (LPMC)** dataset, based on travel survey data collected in London between 2012 and 2015.

The dataset contains:

- **81,086 trips**
- **31,954 individuals**
- **32 original variables**
- **Four travel modes:** Walking, Cycling, Public Transport, and Car

It includes observed travel choices, estimated travel times and costs for alternative modes, trip characteristics, and individual and household attributes.

The original data are stored locally in `data/raw/`, while converted and processed files are stored in `data/processed/`. Dataset files are excluded from the GitHub repository.

## Project Structure

```text
travel-choice-modeling/
├── data/
│   ├── raw/
│   └── processed/
├── notebooks/
│   ├── 01_data_exploration.ipynb
│   └── 02_mnl_model.ipynb
├── src/
├── results/
├── README.md
├── requirements.txt
└── .gitignore
```

## Current Progress

### 01 — Exploratory Data Analysis (Completed)

The first notebook explores the structure of the LPMC dataset and examines observed travel behaviour
