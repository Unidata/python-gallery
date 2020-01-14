"""
GUI for Storm Prediction Center
===============================

This GUI takes in SPC data and plots
all events from a day. It is meant
to show the ease of data access and
creation of a useful plot using the
Siphon Simple Web Service for the SPC.
"""


import cartopy.crs as ccrs
import cartopy.feature as cfeature
from IPython.display import display
import ipywidgets as widgets
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
from siphon.simplewebservice import spc


class SPC_GUI:
    """
    Graphic User Interface designed to allow users to access Storm Prediction Center data.

    This class uses ipython widgets. NOTE: date chosen must be from 2011 to present because
    methods of data parsing are different and will not work for dates earlier.
    """
    def __init__(self):
        """
        Create object that references SPC.py and thus the Storm Prediction Center.

        This initiation creates the SPC object and also creates a widget that allows
        the user to select a date for which to animate SPC reports from the years of
        2011 and onward.
        """
        self.datepicker = widgets.DatePicker(description='Pick a Date:', disabled=False)
        widgets.interact(self.format_date, datepicker=self.datepicker)

    def format_date(self, datepicker):
        """
        Allow user to chose a date for which to plot the SPC events.

        Parameters
        ==========
        datepicker: ipywidget
            allows chosing of date
        """
        if datepicker is not None:
            year_string = datepicker.strftime('%Y')
            month_string = datepicker.strftime('%m')
            day_string = datepicker.strftime('%d')
            self.date_string = year_string + month_string + day_string
            self.stormpicker = widgets.SelectMultiple(options=['tornado', 'hail', 'wind'],
                                                      description='Event Type: ')
            widgets.interact(self.fetch_spc, stormpicker=self.stormpicker)

    def fetch_spc(self, stormpicker):
        """
        Use a date chosen by the user and create widgets that allow you to iterate through
        the pandas dataframe holding the data.

        Parameters
        ==========
        stormpicker: ipywidget
            allows choice of the storm type to be plotted
        """
        self.spc_data_table = []
        for storm_type in stormpicker:
            one_storm_type = spc.SPC(storm_type, self.date_string)
            # Storm must be after year 2011
            if int(self.date_string[:4]) < 2011:
                raise ValueError('SPC GUI does not support events before 2012.')
            self.spc_data_table.append(one_storm_type.day_table)

        if self.spc_data_table != []:
            self.plot_button = widgets.ToggleButton(value=False, description='Plot SPC Events',
                                                    disabled=False, button_style='danger',
                                                    tooltip='Description')

            self.all_events = pd.concat(self.spc_data_table, sort=False)
            self.all_events = self.all_events.sort_index()

            if self.all_events.empty:
                raise ValueError('No storm reports of any type for this date.')

            self.plot_slider = widgets.IntSlider(min=0, max=(len(self.all_events)-1),
                                                 value=0, description='Event #',
                                                 disabled=False)
            self.play = widgets.Play(interval=900, min=0, max=(len(self.all_events)-1),
                                     value=0, description='Event #')
            widgets.jslink((self.plot_slider, 'value'), (self.play, 'value'))
            self.plot_button = widgets.ToggleButton(value=False,
                                                    description='Plot SPC Events',
                                                    disabled=False, button_style='danger',
                                                    tooltip='description')
            self.box = widgets.HBox([self.plot_slider, self.play])
            display(self.plot_slider)
            widgets.interact(self.plot_events, plot_slider=self.play,
                             plot_button=self.plot_button)

    def plot_events(self, plot_slider, plot_button):
        """
        Plot the storm events chosen for the date picked by the user and the type of storms
        selected.

        Parameters
        ==========
        plot_sider: ipywidget
            pases a iteration of the dataframe holding all storm data to be plotted

        plot_button: ipywidget
            toggle button that when activated allows for plotting
        """
        spc_event = self.all_events.iloc[plot_slider]
        if self.plot_button.value is True:
            trimmed_event = spc_event.dropna()
            # Plotting the tracks on top of a cartopy stock image projection
            self.fig = plt.figure(figsize=(14, 11))
            self.ax = self.fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
            self.ax.stock_img()
            state_lines = 'admin_1_states_provinces_lines'
            states_provinces = cfeature.NaturalEarthFeature(category='cultural',
                                                            name=state_lines,
                                                            scale='50m',
                                                            facecolor='none')
            self.ax.add_feature(cfeature.BORDERS)
            self.ax.add_feature(cfeature.COASTLINE)
            self.ax.add_feature(cfeature.LAKES)

            self.ax.add_feature(states_provinces, edgecolor='gray')
            self.ax.set_title('SPC Events for {}/{}/{}'.format(self.date_string[4:6],
                                                               self.date_string[6:8],
                                                               self.date_string[0:4]))
            self.data_projection = ccrs.PlateCarree()
            if 'Size (in)' in trimmed_event:
                self.ax.plot(spc_event['Lon'], spc_event['Lat'], marker='o', color='blue',
                             label='Hail', transform=self.data_projection, markersize=7)
            if 'Speed (kt)' in trimmed_event:
                self.ax.plot(spc_event['Lon'], spc_event['Lat'], marker='o', color='green',
                             label='Wind', transform=self.data_projection, markersize=7)
            if 'F-Scale' in trimmed_event:
                self.ax.plot(spc_event['Lon'], spc_event['Lat'], marker='o', color='red',
                             label='Tornado', transform=self.data_projection, markersize=7)
            self.ax.set_extent([-60, -130, 23, 50])
            hail_event = mpatches.Patch(color='blue', label='Hail')
            wind_event = mpatches.Patch(color='green', label='Wind')
            torn_event = mpatches.Patch(color='red', label='Tornado')
            self.ax.legend(handles=[hail_event, wind_event, torn_event])
            comment = spc_event['Comment']
            if len(comment) > 70:
                comment = comment[:70] + '\n' + comment[70:]
                if len(comment) > 140:
                    comment = comment[:140] + '\n' + comment[140:]
            self.fig.text(0.15, 0.2, comment, fontsize=15)


######################################################################
SPC_GUI()
