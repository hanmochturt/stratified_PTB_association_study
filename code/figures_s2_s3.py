#!/usr/bin/python3
"""
Hannah Takasuka, rotation in the Capra Lab
January 2023
"""

import analysis
import demographics_plots as dp
import conditions_dates_analysis_plots as cd
import stratify_analysis


import matplotlib
import matplotlib.pyplot as plt


def main() -> None:
    matplotlib.use('Qt5Agg')
    matplotlib.rcParams['font.sans-serif'] = "Arial"  # set default font
    cohort = dp.Cohort()
    cohort_cd = cd.get_cohort_patients_conds_dates()
    fig = plt.figure()
    for i, label in enumerate(('A', 'B', 'C')):
        ax = fig.add_subplot(2, 2, i + 1)
        ax.text(-0.1, 1.1, label, transform=ax.transAxes,
                fontsize=16, fontweight='bold', va='top', ha='right')
        if label == 'A':
            cohort.plot_spont_ind_box_whisker_violin(analysis.AGE, ax=ax, title="Maternal Age",
                                                     ylabel="maternal age (years)")
        elif label == 'B':
            cohort.plot_spont_ind_box_whisker_violin(dp.PHECODES_PER_PATIENT, ax=ax, outliers=False,
                                                     title="Count of Phenotypes per Patient",
                                                     ylabel="phenotypes per patient (n)")
        elif label == 'C':
            stratified_cohort = stratify_analysis.stratify_spont_ind_preterm(cohort_cd, term_separate=True)
            variable = cd.FIRST_VISIT
            data = {"indicated": list(stratified_cohort["indicated"][variable]),
                    "spontaneous": list(stratified_cohort["spontaneous"][variable]),
                    "term": list(stratified_cohort["term"][variable])}
            dp.plot_box_whisker_violin_ax(data, ax, title="Time in EHR before Conception",
                                          ylabel="years before conception", outliers=False)

    fig2 = plt.figure()
    for i, label in enumerate(('A', 'B')):
        ax = fig2.add_subplot(2, 1, i + 1)
        ax.text(-0.1, 1.1, label, transform=ax.transAxes,
                fontsize=16, fontweight='bold', va='top', ha='right')
        if label == 'A':
            dp.plot_diagnosis_frequency_histogram(cohort.cohort, ax=ax)
        elif label == 'B':
            cd.plot_diagnoses_time_before_conception_stratified_hist(cohort_cd, ax=ax, log=False)
    plt.show()


if __name__ == "__main__":
    main()
