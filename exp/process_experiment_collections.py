# -*- coding: utf-8 -*-
# Copyright (c) 2021 Wageningen Environmental Research, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), March 2021
"""Processes the YAML files that define the experiment collections in the 'exp' folder

Each experiment is run and the summary results are collected for each experiment in order
to summarize the results over the entire experiment collection.
"""
from pathlib import Path
import importlib
import math
import pickle
import tempfile

import yaml
import pandas as pd
import matplotlib
pd.plotting.register_matplotlib_converters()
matplotlib.use("Agg")
import matplotlib.style
matplotlib.style.use("ggplot")
import matplotlib.pyplot as plt

from pcse.base import ParameterProvider, WeatherDataProvider, WeatherDataContainer

this_dir = Path(__file__).parent
exp_dir = this_dir
exp_results_dir = None


class YAMLWeatherDataProvider(WeatherDataProvider):
    """A WeatherDataProvider which loads the weatherdata contained within the YAML file
        """

    def __init__(self, yaml_weather):
        super().__init__()
        for weather in yaml_weather:
            if "SNOWDEPTH" in weather:
                weather.pop("SNOWDEPTH")
            wdc = WeatherDataContainer(**weather)
            self._store_WeatherDataContainer(wdc, wdc.DAY)


def run_model(yaml_exp):
    """Runs the model for given experiment and returns the simulation results.
    """
    # Load weather, parameters and agromanagement
    wdp = YAMLWeatherDataProvider(yaml_exp["WeatherVariables"])
    params = ParameterProvider(cropdata=yaml_exp["ModelParameters"])
    agro = yaml_exp["Agromanagement"]

    # Retrieve the model we need to run
    module = importlib.import_module(yaml_exp["Model"]["module"])
    model_class = getattr(module, yaml_exp["Model"]["model"])

    # Start and run the model
    model = model_class(params, wdp, agro)
    model.run_till_terminate()

    # Retrieve and return output
    df = pd.DataFrame(model.get_output())
    df = df.set_index("day")
    df_summary = pd.DataFrame(model.get_summary_output())

    return df, df_summary


def get_experimental_data(yaml_exp):
    """gets the experimental data from the YAML file and returns the time-series and summary variables as a dataframe.
    """
    df_exp = pd.DataFrame(yaml_exp["TimeSeriesObservations"]).set_index("day")
    df_exp_sum = pd.DataFrame(yaml_exp["SummaryObservations"])
    return df_exp, df_exp_sum


def build_experiment_plots(df, df_exp, metadata):
    """Plots the experimental observed variables against the simulated variables
    """

    # Determine how many plots we need, depending on number of variables in the experiment
    nfigs = len(df_exp.columns)
    if nfigs == 1:
        fig, axes = plt.subplots(figsize=(12,8))
        axes = [axes,]
    else:
        nrows = math.ceil(nfigs/2.)
        fig, axes = plt.subplots(nrows=nrows, ncols=2, figsize=(12,6*nrows))
        axes = axes.flatten()

    # Loop over variables and axes to make plots
    for varname, ax in zip(df_exp.columns, axes):
        variable_info = metadata["variables"][varname]
        ax.plot(df_exp[varname], marker="o", color="red", linestyle="", label="Observed")
        ax.plot(df[varname], marker="", linestyle="-", color="black", label="Simulated")
        ax.set_title(variable_info["short_title"])
        ax.set_ylabel(variable_info["unit"])
        ax.legend()
    title = "{location} - {year}: {experiment_type}".format(**metadata)
    fig.suptitle(title)

    fname = "{crop}_{location}_{year}_{experiment_type}.png".format(**metadata)
    fname_fp = exp_results_dir / fname
    fig.savefig(fname_fp)
    plt.close("all")


def load_experiment(fname):
    """Load the experiment file and create a cache file if needed
    """
    f = exp_dir / fname
    f_cache = Path(str(f) + ".cache")
    if f_cache.exists():
        if f_cache.stat().st_mtime > f.stat().st_mtime:
            with open(f_cache, "rb") as fp:
                r = pickle.load(fp)
            return r

    # If no or outdated cache, (re)load YAML file
    r = yaml.safe_load(open(f))
    with open(f_cache, "wb") as fp:
        pickle.dump(r, fp, protocol=pickle.HIGHEST_PROTOCOL)

    return r


def run_experiment(exp_file):
    """Runs the experiment for given experiment file.

    It returns dataframes containing the summary variables for both the experiment and the simulation.
    """
    print(f"  - Processing experiment: {exp_file}")
    yaml_exp = load_experiment(exp_file)
    df, df_summary = run_model(yaml_exp)
    df_exp, df_exp_sum = get_experimental_data(yaml_exp)
    build_experiment_plots(df, df_exp, yaml_exp["Metadata"])
    return df_summary, df_exp_sum


def build_summary_plots(collection_fname, exp_collection, df_summary, df_exp_summary):
    """Builds the plots which summarize the experiment collection.
    """
    nfigs = len(exp_collection["VariablestoSummarize"])
    if nfigs == 1:
        fig, axes = plt.subplots(figsize=(12,8))
        axes = [axes,]
    else:
        nrows = math.ceil(nfigs/2.)
        fig, axes = plt.subplots(nrows=nrows, ncols=2, figsize=(12,6*nrows))
        axes = axes.flatten()

    for varname, ax in zip(exp_collection["VariablestoSummarize"], axes):
        variable_info = exp_collection["Metadata"]["variables"][varname]
        ax.scatter(df_summary[varname], df_exp_summary[varname], marker="o")
        ax.plot(variable_info["vrange"], variable_info["vrange"], marker="", linestyle="--", color="k")
        ax.set_title(f'{variable_info["long_title"]} - [{variable_info["unit"]}]')
        ax.set_ylabel("Observed")
        ax.set_xlabel("Simulated")
        ax.set_xlim(variable_info["vrange"])
        ax.set_ylim(variable_info["vrange"])
        ax.set_aspect("equal")

    exp_basename = collection_fname.name.split(".")[0]
    meta = exp_collection["Metadata"]
    fname = "{base}_{crop}_{experiment_type}.png".format(base=exp_basename, **meta)
    fname_fp = exp_results_dir / fname
    fig.savefig(fname_fp)


def run_experiment_collection(collection_fname):
    """Runs all experiments in given experiment collection.
    """

    exp_collection = yaml.safe_load(open(collection_fname))
    meta = exp_collection["Metadata"]
    print("Processing collection for {crop}: {experiment_type}".format(**meta))

    df_summary = pd.DataFrame()
    df_exp_summary = pd.DataFrame()
    for i, fname in enumerate(exp_collection["Experiments"]):
        df1, df2 = run_experiment(fname)
        df_summary = df_summary.append(df1)
        df_exp_summary = df_exp_summary.append(df2)

    build_summary_plots(collection_fname, exp_collection, df_summary, df_exp_summary)


def main(output_dir=None):

    global exp_results_dir
    if output_dir is None:
        exp_results_dir = Path(tempfile.gettempdir()) / "exp_results"
    else:
        exp_results_dir = Path(output_dir) / "exp_results"
    exp_results_dir.mkdir(parents=True, exist_ok=True)
    print(f"Writing expriment results to: {exp_results_dir}")

    for fname in sorted(exp_dir.iterdir()):
        if fname.suffix != ".yaml":
            continue
        run_experiment_collection(fname)


if __name__ == "__main__":
    main()