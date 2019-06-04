"""
Hurricane Tracker with NHC Data
===============================

By: Aodhan Sweeney

This program is a recreation of the 2014 hur_tracker.py
originally written by Unidata Intern Florita Rodriguez. The
2019 version comes with updated interface and functionality,
as well as changing certain dependencies.

"""

import gzip
from io import BytesIO
from io import StringIO

import cartopy.crs as ccrs
import ipywidgets as widgets
import matplotlib.pyplot as plt
import numpy as np
from pandas import DataFrame
import requests


def readUrlFile(url):
    """readUrlFile is a function created to read a .dat file from a given url
    and compile it into a list of strings with given headers."""

    headers = {'User-agent': 'Unidata Python Client Test'}
    response = requests.get(url, headers=headers)
    # Store data response in a string buffer
    string_buffer = StringIO(response.text)
    # Read from the string buffer as if it were a physical file
    data = string_buffer.getvalue()
    return data.splitlines()


def readGZFile(url):
    """readGZFile is a function which opens and reads zipped files. In this case it takes in a
    .gzfile containing information on each storm and returns a byte buffer split based on
    lines."""

    headers = {'User-agent': 'Unidata Python Client Test'}
    response = requests.get(url, headers=headers)
    # Store data response in a bytes buffer
    bio_buffer = BytesIO(response.content)
    # Read from the string buffer as if it were a physical file
    gzf = gzip.GzipFile(fileobj=bio_buffer)
    data = gzf.read()
    return data.splitlines()


def split_storm_info(storm_list):
    """split_storm_info takes a list of strings and creates a pandas dataframe
    for the data set taken off the NHC archive. This function is called in the main to
    find all storms."""

    name, cycloneNum, year, stormType, basin, filename = [], [], [], [], [], []
    for line in storm_list[1:]:
        fields = line.split(',')
        name.append(fields[0].strip())
        basin.append(fields[1].strip())
        cycloneNum.append(fields[7].strip())
        year.append(fields[8].strip())
        stormType.append(fields[9].strip())
        filename.append(fields[-1].strip().lower())

    storms = DataFrame({'Name': name, 'Basin': basin, 'CycloneNum': np.array(cycloneNum),
                        'Year': np.array(year), 'StormType': stormType,
                        'Filename': filename})
    return(storms)


class Storm_Selection_gui:
    """Storm_Selection_gui is a graphic user interface object designed for selecting storms
    from the National Hurricane Center's storm list database."""

    def __init__(self):
        """__init__ is the initiation function for the Storm_Selection_gui object. This
        initiation creates a dataframe for the storm data from the National Hurricane
        Center. In addition, this initiation creates widgets to select a specific storm
        from the National Hurricane Center database with both a widget for a year
        selection and also for a track button which actually retrieves the models
        and tracks for a given storm. """

        # Setting up storm object table
        fileLines = readUrlFile('http://ftp.nhc.noaa.gov/atcf/index/storm_list.txt')
        self.storm_table = split_storm_info(fileLines)

        # Creation of widgets
        # Year Selector Slider (year_slider)
        self.year_slider = widgets.IntSlider(min=1851, max=2019, value=2019,
                                             description='Storm Year: ')
        widgets.interact(self.get_storms_slider, year_slider=self.year_slider)
        # Button to retrieve storm tracks (track_button)
        self.track_button = widgets.ToggleButton(value=False, description='Get Storm Tracks',
                                                 disabled=False, button_style='info',
                                                 tooltip='Description')
        widgets.interact(self.get_tracks, track_button=self.track_button)

    def get_storms_slider(self, year_slider):
        """get_storms_slider is a function that is written to take a user defined year from
        1851 to 2019 and then trim the split_storms_info dataframe to just have
        storms from that year. This takes in the interactive widget called year_slider."""

        self.one_year_table = self.storm_table[self.storm_table.Year == str(year_slider)]
        self.storm_names = widgets.Dropdown(options=self.one_year_table['Name'],
                                            description='Storm Names: ')
        widgets.interact(self.get_name_dropdown, storm_names=self.storm_names)

    def get_name_dropdown(self, storm_names):
        """get_name_dropdown is a function that allows for selection of a specific storm
        based on its name. This function returns the filename of a given storm. This
        function interacts with the widget called forecast_select."""

        name = self.storm_names.value
        one_storm_row = self.one_year_table[self.one_year_table.Name == name]
        self.filename = one_storm_row.Filename
        file_name = self.filename.tolist()
        if self.filename.empty is False:
            self.filename = file_name[0]
        elif self.filename.empty is True:
            print('No file name data for this year.')

    def get_models(self, models):
        """get_models is a function that is linked to the selection widget which allows for
        multiple models to be selected to plot for a specific hurricane track."""

        self.model_selection_latlon(models)

    def time_slider(self, time_slide):
        """time_slider is a function which changes the time for which to plot the models for
        a given hurricane track."""

        time = self.date_times[time_slide]
        return(time)

    def get_tracks(self, track_button):
        """get_tracks is a function that will create the url and pull the data for either
        the forecast track or best track for a given storm. The Url is made by using both
        the year and the filename. This function will then read the data and create a data
        frame for both the forecast and best tracks and compile these data frames into a
        dictionary. This function returns this dictionary of forecast and best tracks. """

        year = str(self.year_slider.value)
        filename = self.filename
        data_dictionary = {}
        # Current year data is stored in a different location
        if year == '2019':
            urlf = 'http://ftp.nhc.noaa.gov/atcf/aid_public/a{}.dat.gz'.format(filename)
            urlb = 'http://ftp.nhc.noaa.gov/atcf/btk/b{}.dat'.format(filename)
        else:
            urlf = 'http://ftp.nhc.noaa.gov/atcf/archive/{}/a{}.dat.gz'.format(year, filename)
            urlb = 'http://ftp.nhc.noaa.gov/atcf/archive/{}/b{}.dat.gz'.format(year, filename)

        url_links = [urlf, urlb]
        url_count = 0
        if track_button is True:
            for url in url_links:
                # Checking if url is valid, if status_code is 200 then website is active
                if requests.get(url).status_code == 200:
                    if url.endswith('.dat'):
                        lines = readUrlFile(url)
                    else:
                        lines = readGZFile(url)

                    # Splitting the method for which we will create the dataframe
                    lat, lon, basin, cycloneNum = [], [], [], []
                    warnDT, model, forecast_hour = [], [], []
                    for line in lines:
                        line = str(line)
                        line = line[2:]
                        fields = line.split(',')
                        # Joins together lattitude and longitude strings without
                        # directional letters.
                        # Includes type conversion in order to divide by 10 to
                        # get the correct coordinate.
                        latSingle = int(fields[6][:-1])/10.0
                        lonSingle = -(int(fields[7][:-1])/10.0)
                        lat.append(latSingle)
                        lon.append(lonSingle)
                        basin.append(fields[0])
                        forecast_hour.append(fields[5])
                        cycloneNum.append(fields[1].strip())
                        warnDT.append(fields[2].strip())
                        model.append(fields[4].strip())

                        # Combining data from file into a Pandas Dataframe.
                        storm_data_frame = DataFrame({'Basin': basin,
                                                      'CycloneNum': np.array(cycloneNum),
                                                      'WarnDT': np.array(warnDT),
                                                      'Model': model, 'Lat': np.array(lat),
                                                      'Lon': np.array(lon),
                                                      'forecast_hour':
                                                      np.array(forecast_hour)})
                        # Adding this newly created DataFrame to a dictionary
                        if url_count == 0:
                            data_dictionary['forecast'] = storm_data_frame
                        else:
                            data_dictionary['best_track'] = storm_data_frame

                else:
                    print('url {} was not valid, select different storm.'.format(url))
                    track_button = False

                url_count += 1

            self.storm_dictionary = data_dictionary
            forecast = data_dictionary.get('forecast')

            unique_models, unique_index = list(np.unique(forecast['Model'].values,
                                               return_index=True))

            # Selection tool to pick from available models to plot (model_select)
            self.model_select = widgets.SelectMultiple(options=unique_models,
                                                       value=[unique_models[0]],
                                                       description='Models: ',
                                                       disabled=False)
            widgets.interact(self.get_models, models=self.model_select)

    def model_selection_latlon(self, models):
        """model_selection_latlon is a function that allows the user to select a model for
        a given storm and whether the tracks are forecast or best tracks. The parameters
        for this are a string stating whether the user wants forecast or best tracks and
        also all model outputs for all forecasts and best tracks compiled into a python
        dictionary. The latlon part of this function comes from taking the users selected
        model and getting the latitudes and longitudes of all positions of the storm for
        this forecast. This function then returns these lats and lons as a pandas.Series"""

        if self.model_select.disabled is False:
            # We will always plot best track, and thus must save the coordinates for plotting
            best_track = self.storm_dictionary.get('best_track')
            self.date_times = best_track['WarnDT']
            lats = best_track['Lat']
            lons = best_track['Lon']
            self.best_track_coordinates = [lats, lons]

            model_tracks = self.storm_dictionary.get('forecast')

            self.model_table = []
            for model in models:
                one_model_table = model_tracks[model_tracks['Model'] == model]
                self.model_table.append(one_model_table)

            # Slider to choose time frame for which to plot models (plot_slider)
            self.plot_slider = widgets.IntSlider(min=0, max=(len(self.date_times)-1),
                                                 value=0, description='Tracks Time',
                                                 disabled=False)
            widgets.interact(self.plotting, plot_slider=self.plot_slider)

    def plotting(self, plot_slider):
        """plotting is a function that that plots the models tracks and best tracks for all
        selected models. These tracks are then plotted on the Blue Marble projection. """

        if self.plot_slider.disabled is False:

            # Identifying the time associated with the models for time text box
            year = self.date_times[plot_slider][0: 4]
            month = self.date_times[plot_slider][4: 6]
            day = self.date_times[plot_slider][6: 8]
            hour = self.date_times[plot_slider][8: 10]
            time_string = 'Date: {0}/{1}/{2} \n Hour: {3}'.format(month, day, year, hour)

            # Finding data for best track, and extremes for which to base axis extent on
            self.best_lats = np.array(self.best_track_coordinates[0])
            self.best_lons = np.array(self.best_track_coordinates[1])
            min_best_lat = min(self.best_lats)
            max_best_lat = max(self.best_lats)
            min_best_lon = min(self.best_lons)
            max_best_lon = max(self.best_lons)

            # Plotting the track on a cartopy stock image
            self.fig = plt.figure(figsize=(14, 11))
            self.ax = self.fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
            self.ax.stock_img()

            self.data_projection = ccrs.PlateCarree()
            self.ax.plot(self.best_lons, self.best_lats, marker='o', color='white',
                         label='Best Track', transform=self.data_projection)
            self.ax.set_extent([(min_best_lon - 30), (max_best_lon + 30),
                                (min_best_lat - 30), (max_best_lat + 30)])

            jet = plt.get_cmap('jet')
            colors = iter(jet(np.linspace(0.2, 1, (len(self.model_select.value)+1))))
            left = .1
            bottom = .1
            self.ax.text(left, bottom, time_string, transform=self.ax.transAxes,
                         fontsize=14, color='black')

            for model_type in self.model_table:
                one_model_time = model_type[model_type['WarnDT'] ==
                                            self.date_times[plot_slider]]
                lats = one_model_time['Lat'].tolist()
                lons = one_model_time['Lon'].tolist()
                if len(lats) != 0:
                    model_list = model_type['Model'].tolist()
                    self.ax.plot(lons, lats, marker='o', color=next(colors),
                                 label=model_list[0])

            plt.title('Storm Name: {0} Year: {1}'.format(self.storm_names.value,
                      str(self.year_slider.value)))
            plt.legend()


######################################################################
storm_selection = Storm_Selection_gui()
