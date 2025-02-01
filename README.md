# stratified_PTB_association_study

I copied Jackie's RPL readme; work in progress that I am updating

This repository contains the code associated with our paper:

Pre-conception clinical risk factors differ between spontaneous and indicated preterm birth in a densely phenotyped EHR cohort

Jean M. Costello†, Hannah Takasuka†, Jacquelyn Roger, Ophelia Yin, Alice Tang, Tomiko Oskotsky, Marina Sirota*, and John A. Capra*
†co-first authors
*co-senior authors

Part of our study utilized electronic health record (EHR) data that was standardized to the Observational Medical Outcomes Partnership (OMOP) common data model. To replicate this study using another OMOP EHR database, follow the query instructions below and then use the code in these R markdown scripts to filter your patients, aggregate diagnoses, implement analyses, analyze results, and create figures.

The other part of our study utilized data column names that are specific to the University of California, San Francisco (i.e., clinician-curated perinatal database (PDB) and Clinical Data Warehouse EHR).

Query instructions:
1. Query your OMOP database to identify patients that meet the initial inclusion criteria for the RPL group. The initial RPL criteria is: any record of pregnancy loss (as defined in the pregnancy loss, recurrent pregnancy loss, and history of pregnancy loss sections of Supplementary File 1).
2. Query your OMOP database to identify patients that meet the initial inclusion criteria for the Control group. The initial Control criteria is: any record of live-birth (as defined in the live-birth section of Supplementary File 1).
3. Query the included RPL and Control patients' subsets of the following OMOP tables: person, condition_occurrence, observation, procedure_occurrence, and visit_occurrence.
4. The inputs for Rmd 01 are: the table subsets from #3, the OMOP concept table, and the concept lists in Supplementary File 1.

System requirements: This code was developed and tested using R version 4.0.2 on the macOS Monterey version 12.6.3.

Installation guide: To run the code, clone this repository. Cloning typically takes a minute or less.

Code and Data Sharing Permissions

The University of California, San Francisco electronic health record data used for our project is not publicly available. If you would like to set up an official collaboration, email marina.sirota@ucsf.edu. If you are a UCSF affiliate with OMOP and CDW and/or PDB access, email hannah.takasuka@ucsf.edu to request the version of this repository with EHR data and exact CDW and PDB queries. Requests should be processed within a couple of weeks.

License:

    Our RPL association study code analyzes EHR data from RPL and control patients, and identifies diagnoses associated with RPL.
    Copyright (C) 2023 Jacquelyn Roger

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: jacquelyn.roger@ucsf.edu
    Located at Bakar Institute at UCSF
