"""
GUI for National Hurricane Center Data
======================================
By: Aodhan Sweeney

This GUI takes data pulled from the National
Hurricane Center and plots specific storms and
their tracks. The GUI runs using ipywidgets, and
is intended to be used as an Jupyter notebook.
"""

import cartopy.crs as ccrs
from IPython.display import display
import ipywidgets as widgets
import matplotlib.pyplot as plt
import numpy as np
from siphon.simplewebservice import nhc


class NHC_GUI:
    """
    Graphic User Interface designed to allow users to access National Hurricane Center data.

    This class uses ipython widgets, and the order in which the functions appear in this script
    correspond to the order in which the functions and widgets are called/used.
    """

    def __init__(self):
        """
        Create object that references NHC.py and thus the National Hurricane Center.

        This initiation creates the National Hurricane Center object and also creates
        a widget that allows the user to 1.) select the year in which to find storms,
        and 2.) indicate when they have chosen the year and when to continue with parsing
        the storm_table.
        """
        self.NHCD = nhc.NHCD()
        self.storm_table = self.NHCD.storm_table
        # Year Slider Widget to select year for which to retrieve storm data.
        self.year_slider = widgets.IntSlider(min=1851, max=2019, value=2019,
                                             description='Storm Year: ')
        widgets.interact(self.get_storms_slider, year_slider=self.year_slider)
        # Storm Track toggle button to initiate storm track retrieval.
        self.track_button = widgets.ToggleButton(value=False, description='Get Storm Tracks',
                                                 disabled=False, button_style='info',
                                                 tooltip='Description')
        widgets.interact(self.get_track, track_button=self.track_button)

    def get_storms_slider(self, year_slider):
        """
        Take a year chosen by the user, and create a list of all storms during that year.

        Parameters
        ----------
        year_slider: ipywidget
            tells value of year chosen by the user
        """
        self.year = str(year_slider)
        self.one_year_table = self.storm_table[self.storm_table.Year == year_slider]
        self.storm_names = widgets.Dropdown(options=self.one_year_table['Name'],
                                            description='Storm Names: ')
        widgets.interact(self.get_name_dropdown, storm_names=self.storm_names)

    def get_name_dropdown(self, storm_names):
        """
        Take names created by previous function and allow selection of which models to plot.

        Parameters
        ----------
        storm_names: list
            all storms in a given year
        """
        name = self.storm_names.value
        one_storm_row = self.one_year_table[self.one_year_table.Name == name]
        self.filename = one_storm_row.Filename
        file_name = self.filename.tolist()
        if self.filename.empty is False:
            self.filename = file_name[0]
        elif self.filename.empty is True:
            raise ValueError('No data for specific file name of the chosen.')

    def get_track(self, track_button):
        """
        Query whether track button has been toggled, and create select model widget.

        Parameters
        ----------
        track_button: ipywidget
            button that when toggled indicates that user is ready to select the model
            tracks for the chosen storm.
        """
        if self.track_button.value is True:
            unique_models = self.NHCD.get_tracks(self.year, self.filename)
            self.model_select = widgets.SelectMultiple(options=unique_models,
                                                       value=[unique_models[0]],
                                                       description='Models: ',
                                                       disabled=False)
            widgets.interact(self.get_models, models=self.model_select)

    def get_models(self, models):
        """
        Select models from a list of all models for a given stormself.

        Parameters
        ----------
        models: list
            list of ran for a given storm
        """
        self.NHCD.model_selection_latlon(models)
        self.date_times = self.NHCD.date_times.tolist()
        self.plot_slider = widgets.IntSlider(min=0, max=(len(self.date_times)-1),
                                             value=0, description='Tracks Time',
                                             disabled=False)
        self.play = widgets.Play(interval=800, min=0, max=(len(self.date_times)-1),
                                 value=0, description='Tracks Time')
        widgets.jslink((self.plot_slider, 'value'), (self.play, 'value'))
        self.box = widgets.HBox([self.plot_slider, self.play])
        display(self.plot_slider)
        widgets.interact(self.plotting, plot_slider=self.play)

    def plotting(self, plot_slider):
        """
        Plot selected model tracks and best track for given storm.

        Parameters
        ----------
        plot_slider: ipywidget
            Widget where assigned value is plot/play widget value
        """
        if self.plot_slider.disabled is False:

            # Identifying the time associated with the models for time text box
            year = self.date_times[plot_slider][0: 4]
            month = self.date_times[plot_slider][4: 6]
            day = self.date_times[plot_slider][6: 8]
            hour = self.date_times[plot_slider][8: 10]
            time_string = 'Date: {0}/{1}/{2} \nHour: {3}'.format(month, day, year, hour)

            # Finding data for best track, and extremes for which to base axis extent on
            self.best_lats = np.array(self.NHCD.best_track_coordinates[0])
            self.best_lons = np.array(self.NHCD.best_track_coordinates[1])
            min_best_lat = min(self.best_lats)
            max_best_lat = max(self.best_lats)
            min_best_lon = min(self.best_lons)
            max_best_lon = max(self.best_lons)

            # Plotting the tracks on top of a cartopy stock image projection
            self.fig = plt.figure(figsize=(14, 11))
            self.ax = self.fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
            self.ax.stock_img()
            self.data_projection = ccrs.PlateCarree()
            self.ax.plot(self.best_lons, self.best_lats, marker='o', color='white',
                         label='Best Track', transform=self.data_projection)

            self.ax.set_extent([(min_best_lon - 40), (max_best_lon + 40),
                                (min_best_lat - 40), (max_best_lat + 40)])

            jet = plt.get_cmap('jet')
            colors = iter(jet(np.linspace(0.2, 1, (len(self.model_select.value)+1))))
            left = .1
            bottom = .1
            self.ax.text(left, bottom, time_string, transform=self.ax.transAxes,
                         fontsize=14, color='black')

            for model_type in self.NHCD.model_table:
                one_time = model_type[model_type['WarnDT'] == self.date_times[plot_slider]]
                lats = one_time['Lat'].tolist()
                lons = one_time['Lon'].tolist()
                if len(lats) != 0:
                    model_list = model_type['Model'].tolist()
                    self.ax.plot(lons, lats, marker='o', color=next(colors),
                                 label=model_list[0])
            plt.title('Storm Name: {0} Year: {1}'.format(self.storm_names.value,
                                                         self.year))
            plt.legend()


######################################################################
NHC_GUI()
