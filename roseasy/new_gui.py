#!/usr/bin/env python2
# encoding: utf-8

"""\
Judge forward-folded candidates in computational protein design pipelines.

Usage:
    show_my_designs.py [options] <pdb_directories>...
    show_my_designs.py --version

Options:
    -F, --no-fork
        Do not fork into a background process.

    -f, --force
        Force the cache to be regenerated.

    -q, --quiet
        Build the cache, but don't launch the GUI.

    -v, --version
        Print the version number and exit.

Features:
    1. Extract quality metrics from forward-folded models and plot them against
       each other in any combination.

    2. Easily visualize specific models by right-clicking on plotted points.
       Add your own visualizations by writing `*.sho' scripts.

    3. Plot multiple designs at once, for comparison purposes.

    4. Keep notes on each design, and search your notes to find the designs you
       want to visualize.
"""

## Imports
from gi import pygtkcompat

pygtkcompat.enable()
pygtkcompat.enable_gtk(version='3.0')
import collections, glob, gzip, os, re, shutil, subprocess, sys
import matplotlib
matplotlib.use('GTK3Agg')
from gi.repository import Gtk as gtk
from gi.repository import Pango as pango
from gi.repository import GObject as gobject
from gi.repository import Gdk as gdk
import yaml
import matplotlib.pyplot as plt, numpy as np, scipy as sp, pandas as pd, seaborn as sns

from matplotlib.figure import Figure
from matplotlib.backends.backend_gtk3agg import FigureCanvasGTK3Agg
from matplotlib.backends.backend_gtk3 import NavigationToolbar2GTK3
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
from pprint import pprint
from roseasy import structures


class Design (object):

    def __init__(self, directory, use_cache=True):
        self.directory = directory
        self.cache_path = os.path.join(directory, 'models.pkl')
        self.notes_path = os.path.join(directory, 'notes.txt')
        self.rep_path = os.path.join(directory, 'representative.txt')

        self._models = None
        self._metrics = {}
        self._notes = ""
        self._representative = None

        self._load_models(use_cache)
        self._load_annotations()

    def __str__(self):
        return '<ModelGroup dir={}>'.format(self.directory)

    def __len__(self):
        return len(self.paths)


    @property
    def paths(self):
        return self._models['path']

    @property
    def notes(self):
        return self._notes

    @notes.setter
    def notes(self, notes):
        self._notes = notes
        self._save_notes()

    @property
    def representative(self):
        if self._representative is None:
            return self.get_metric('total_score').idxmin()
        else:
            return self._representative

    @representative.setter
    def representative(self, index):
        self._representative = index
        self._save_representative()

    @property
    def representative_path(self):
        return self.paths[self.representative]

    @property
    def metrics(self):
        return self._metrics

    def get_metric(self, metric):
        if metric not in self.metrics:
            message = "No such metric: '{}'\n".format(metric)
            message += "Defined metrics are: " + ', '.join(
                    "'{}'".format(x) for x in self.metrics)
            print(type(metric), ' '.join(str(type(x)) for x in self.metrics))

            raise RuntimeError(message)

        return self._models[metric]

    def get_coord(self, x_metric, y_metric, index=None):
        i = index if index is not None else self.representative
        return self.get_metric(x_metric)[i], self.get_metric(y_metric)[i]


    def _load_models(self, use_cache):
        """
        Load a variety of score and distance metrics for the structures found
        in the given directory.  As much information as possible will be
        cached.  Note that new information will only be calculated for file
        names that haven't been seen before.  If a file changes or is deleted,
        the cache will not be updated to reflect this and you may be presented
        with stale data.
        """

        self._models, self._metrics = structures.load(
                self.directory,
                use_cache=use_cache,
                require_io_dir=False,
        )

    def _load_metrics(self):
        # Treat column in self._models that contains numeric data as a metric.
        # Any dtype other than 'object' is assumed to be numeric.

        self._metrics = {
            x: MetricInfo(
                x,
                title=get_metric_title(x, self),
                order=get_metric_order(x, self),
                guide=get_metric_guide(x, self),
                limits=get_metric_limits(x, self),
            )
            for x in list(self._models.keys())
            if self._models[x].dtype != 'object'
        }

        # Make sure at least two metrics have been associated with each model
        # in this directory.

        if len(self._metrics) == 0:
            raise IOError("no metrics defined for the models in '{}'".format(self.directory))
        if len(self._metrics) == 1:
            name = next(iter(self._metrics))
            raise IOError("only found one metric '{}' for the models in '{}', need at least two".format(name, self.directory))

    def _load_annotations(self):
        try:
            with open(self.notes_path) as file:
                self._notes = file.read()
        except IOError:
            pass

        try:
            with open(self.rep_path) as file:
                self._representative = int(file.read())
        except IOError:
            pass

    def _save_notes(self):
        with open(self.notes_path, 'w') as file:
            file.write(self.notes)

        if os.path.exists(self.notes_path) and not self.notes:
            os.remove(self.notes_path)

    def _save_representative(self):
        if self._representative is not None:
            with open(self.rep_path, 'w') as file:
                file.write(str(self._representative))

        elif os.path.exists(self.rep_path):
            os.remove(self.rep_path)


class ShowMyDesigns (gtk.Window):

    def __init__(self, designs):

        # Setup the parent class.

        gtk.Window.__init__(self)
        self.add_events(gdk.EventMask.KEY_PRESS_MASK)
        self.connect('key-press-event', self.on_hotkey_press)
        self.set_icon_from_file(os.path.join(
            os.path.dirname(os.path.realpath(__file__)), 'icon.png'))

        # Setup the data members.

        self.designs = designs

        self.keys = list()
        self.selected_model = None
        self.is_legend_visible = False
        self.is_representative_visible = False
        self.is_model_count_visible = False

        self.metrics = {
                k: next(iter(self)).metrics[k] \
                for k in set.intersection(*[set(x.metrics) for x in self])
        }
        print(self.metrics)
        self.sorted_metrics = sorted(
                self.metrics,
                key=lambda k: (self.metrics[k].order, self.metrics[k].title)
        )
        print(self.sorted_metrics)
        self.x_metric = (
                default_x_metric
                if default_x_metric in self.metrics
                else self.sorted_metrics[0])
        self.y_metric = (
                default_y_metric
                if default_y_metric in self.metrics
                else self.sorted_metrics[1])

        # Setup the GUI.

        self.connect('destroy', lambda x: gtk.main_quit())
        self.set_default_size(int(1.618 * 630), 630)

        menu_bar = self.setup_menu_bar()
        model_viewer = self.setup_model_viewer()

        vbox = gtk.VBox()
        vbox.pack_start(menu_bar, False)
        vbox.pack_end(model_viewer, True)

        self.add(vbox)
        self.set_border_width(3)
        self.update_everything()
        self.show_all()
        self.set_focus(None)

        n = len(self.designs)

        self.hide_model_list() if n == 1 else self.show_model_list()
        self.hide_filter_pane()
        self.show_annotation_pane()
        self.hide_legend()
        self.show_representative()
        self.hide_model_count()

    def __iter__(self):
        return iter(list(self.designs.values()))


    def setup_menu_bar(self):
        bar = gtk.MenuBar()

        # The "File" menu:
        menu = gtk.Menu()
        item = gtk.MenuItem("_File")
        item.set_submenu(menu)
        bar.append(item)

        item = gtk.MenuItem("Save selected paths")
        item.connect('activate', lambda _: self.save_selected_paths())
        menu.append(item)

        item = gtk.MenuItem("Save selected funnels")
        item.connect('activate', lambda _: self.save_selected_funnels())
        menu.append(item)

        # The "View" menu:
        menu = gtk.Menu()
        item = gtk.MenuItem("_View")
        item.set_submenu(menu)
        bar.append(item)

        item = self.model_list_toggle = gtk.CheckMenuItem("Model list")
        item.connect('activate', self.on_toggle_model_list)
        menu.append(item)

        item = self.filter_pane_toggle = gtk.CheckMenuItem("Filter pane")
        item.connect('activate', self.on_toggle_filter_pane)
        menu.append(item)

        item = self.annotation_pane_toggle = gtk.CheckMenuItem("Annotation pane")
        item.connect('activate', self.on_toggle_annotation_pane)
        menu.append(item)

        item = gtk.SeparatorMenuItem()
        menu.append(item)

        item = self.legend_toggle = gtk.CheckMenuItem("Legend")
        item.connect('activate', self.on_toggle_legend)
        menu.append(item)

        item = self.representative_toggle = gtk.CheckMenuItem("Representative")
        item.connect('activate', self.on_toggle_representative)
        menu.append(item)

        item = self.model_count_toggle = gtk.CheckMenuItem("Model count")
        item.connect('activate', self.on_toggle_model_count)
        menu.append(item)

        return bar

    def setup_model_viewer(self):
        plot = self.setup_plot()
        self.model_list = self.setup_model_list()
        self.filter_pane = self.setup_filter_pane()
        self.annotation_pane = self.setup_annotation_pane()

        self.sidebar = gtk.VPaned()
        self.sidebar.add1(self.model_list)
        self.sidebar.add2(self.filter_pane)

        bottombar = gtk.VPaned()
        bottombar.add1(plot)
        bottombar.add2(self.annotation_pane)

        viewer = gtk.HPaned()
        viewer.add1(self.sidebar)
        viewer.add2(bottombar)

        return viewer

    def setup_model_list(self):
        list_store = gtk.ListStore(str)

        text = gtk.CellRendererText()
        icon = gtk.CellRendererPixbuf()

        self.view = gtk.TreeView(list_store)
        self.view.set_model(list_store)
        self.view.set_rubber_banding(True)
        self.view.set_enable_search(False)
        self.view.set_headers_visible(False)

        columns = [
                ('Name', 'directory'),
        ]

        for index, parameters in enumerate(columns):
            title, attr = parameters

            def cell_data_func(column, cell, model, iter, attr): #
                key = model.get_value(iter, 0)
                design = self.designs[key]
                text = getattr(design, attr)
                cell.set_property('text', text)

            def sort_func(model, iter_1, iter_2, attr): #
                key_1 = model.get_value(iter_1, 0)
                key_2 = model.get_value(iter_2, 0)
                design_1 = self.designs[key_1]
                design_2 = self.designs[key_2]
                value_1 = getattr(design_1, attr)
                value_2 = getattr(design_2, attr)
                return cmp(value_1, value_2)

            list_store.set_sort_func(index, sort_func, attr);

            column = gtk.TreeViewColumn(title, text)
            column.set_cell_data_func(text, cell_data_func, attr)
            column.set_sort_column_id(index)
            self.view.append_column(column)

        selector = self.view.get_selection()
        selector.connect("changed", self.on_select_designs)
        ### If selection is weird, double check this!
        ### gtk.SelectionMode(3) should correspond to GTK_SELECTION_MULTIPLE
        selector.set_mode(gtk.SelectionMode(3))

        scroller = gtk.ScrolledWindow()
        scroller.set_policy(gtk.PolicyType.AUTOMATIC,
                gtk.PolicyType.AUTOMATIC)
        scroller.add(self.view)

        frame = gtk.Frame()
        frame.add(scroller)

        self.search_form = gtk.Entry()
        self.search_form.set_icon_from_stock(gtk.EntryIconPosition.SECONDARY, gtk.STOCK_FIND)

        search_buffer = self.search_form.get_buffer()
        search_buffer.connect('deleted-text', self.on_search_in_notes)
        search_buffer.connect('inserted-text', self.on_search_in_notes)

        vbox = gtk.VBox()
        vbox.pack_start(self.search_form, False, True, 0)
        vbox.pack_start(frame, True, True, 0)

        return vbox

    def setup_plot(self):
        figure = Figure(facecolor='#edecea')

        # Create the axes.

        self.axes = figure.add_axes((0.15, 0.15, 0.75, 0.75))
        self.axes.set_ylabel('Score')

        # Create the canvas.

        self.canvas = FigureCanvas(figure)
        self.canvas.mpl_connect('pick_event', self.on_select_model)
        self.canvas.mpl_connect('button_press_event', self.on_click_plot_mpl)
        self.canvas.mpl_connect('motion_notify_event', self.on_move_mouse_mpl)
        self.canvas.connect('button-press-event', self.on_click_plot_gtk)
        self.canvas.set_size_request(-1, 350)

        # Create the tool bar.

        self.toolbar = NavigationToolbar(self.canvas, self)

        # Place all the widgets.

        self.mouse_position = gtk.Label("")

        table = gtk.Table(3, 5)
        table.attach(self.toolbar, 0, 1, 0, 3)
        table.attach(self.mouse_position, 3, 4, 1, 2, xoptions=0, yoptions=0, xpadding=3)

        vbox = gtk.VBox()
        vbox.pack_start(self.canvas, True, True, 0)
        vbox.pack_start(table, False, True, 0)

        return vbox

    def setup_filter_pane(self):
        pane = FilterPane(self)
        pane.connect('updated', lambda _: self.update_plot())
        return pane

    def setup_annotation_pane(self):
        self.notes = gtk.TextView()
        self.notes.set_wrap_mode(gtk.WRAP_WORD)
        self.notes.set_size_request(-1, 100)
        self.notes.set_left_margin(3)
        self.notes.set_right_margin(3)
        self.notes.set_pixels_above_lines(3)
        self.notes.set_pixels_below_lines(3)
        self.notes.set_cursor_visible(True)
        self.notes.get_buffer().connect('changed', self.on_edit_annotation)

        scroll_window = gtk.ScrolledWindow()
        scroll_window.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)
        scroll_window.add(self.notes)

        frame = gtk.Frame()
        frame.add(scroll_window)

        return frame

    def setup_metric_menu(self, callback=None, initial_choice=None):
        try: self.metric_store
        except AttributeError:
            self.metric_store = gtk.ListStore(str, str)

            for key in self.sorted_metrics:
                metric = self.metrics[key]
                self.metric_store.append([metric.name, metric.title])

        cell = gtk.CellRendererText()
        ### previousoy ComboBox
        menu = gtk.ComboBoxText(model=self.metric_store)
        menu.pack_start(cell, True)
        menu.add_attribute(cell, 'text', 1)

        menu.set_active(0)
        for i, metric in enumerate(self.sorted_metrics):
            if metric == initial_choice:
                menu.set_active(i)

        if callback:
            menu.connect('changed', callback)

        return menu


    def on_hotkey_press(self, widget, event):
        key = gdk.keyval_name(event.keyval).lower()
        if event.state & gdk.CONTROL_MASK: key = 'ctrl-' + key
        if event.state & gdk.SHIFT_MASK: key = 'shift-' + key

        hotkeys = {
                'escape': self.normal_mode,
        }
        normal_mode_hotkeys = {
                'i': self.insert_mode,     'a': self.insert_mode,
                'z': self.zoom_mode,
                'x': self.pan_mode,
                'c': self.refocus_plot,
                'tab': self.cycle_y_metric,
                'space': self.cycle_x_metric,
                'shift-tab': self.reverse_cycle_y_metric,
                'shift-space': self.reverse_cycle_x_metric,
        }
        multi_design_hotkeys = {
                'j': self.next_design,     'f': self.next_design,
                'k': self.previous_design, 'd': self.previous_design,
                'slash': self.search_mode,
        }

        keep_focus = (
                gtk.Entry,
                gtk.TextView,
                gtk.Button,
                gtk.ComboBox
        )
        if not isinstance(self.get_focus(), keep_focus):
            hotkeys.update(normal_mode_hotkeys)
            if len(self.designs) > 1:
                hotkeys.update(multi_design_hotkeys)

        if key in hotkeys:
            hotkeys[key]()
            return True

    def on_search_in_notes(self, entry_buffer, *_):
        self.update_designs()

    def on_select_designs(self, selection):
        new_keys = []
        old_keys = self.keys[:]
        self.keys = []
        model, paths = selection.get_selected_rows()

        for path in paths:
            iter = model.get_iter(path)
            key = model.get_value(iter, 0)
            new_keys.append(key)

        # Don't change the order of designs that were already selected.  The
        # order affects how the color of the design in the score vs rmsd plot,
        # and things get confusing if it changes.

        for key in old_keys:
            if key in new_keys:
                self.keys.append(key)

        for key in new_keys:
            if key not in self.keys:
                self.keys.append(key)

        # Rename the window based on the current selection.

        subtitle = ""
        if len(self.keys) == 1:
            subtitle = " ({})".format(self.keys[0])
        if len(self.keys) > 1:
            subtitle = " ({}, ...)".format(self.keys[0])

        self.set_title("Show My Designs" + subtitle)

        # This is an efficiency thing.  The 'J' and 'K' hotkeys works in two
        # steps: first unselect everything and then select the next row in
        # order.  Redrawing the plot is expensive, so it's worthwhile to skip
        # redrawing after that first step.

        if self.keys:
            self.update_plot()
            self.update_annotations()

    def on_select_model(self, event):
        self.selected_model = event.ind[0], event.artist.design

    def on_move_mouse_mpl(self, event):
        if event.xdata is None or event.ydata is None:
            # The data coordinates will be None only if the mouse is outside
            # the data area.
            self.mouse_position.set_text("")
        else:
            coord = '{:0.2f}, {:0.2f}'.format(event.xdata, event.ydata)
            self.mouse_position.set_text(coord)

    def on_click_plot_mpl(self, event):
        pass

    def on_click_plot_gtk(self, widget, event):
        # Ignore any event that isn't a right button click or a left button
        # click with the control key held down.

        is_right_click = \
                (event.button == 3) or \
                (event.button == 1 and event.get_state() & gdk.CONTROL_MASK)

        if not is_right_click: return
        if self.toolbar._active == 'PAN': return
        if self.toolbar._active == 'ZOOM': return
        if self.selected_model is None: return

        # Figure out which model was clicked.

        index, design = self.selected_model
        rep_index = design.representative
        path = os.path.join(design.directory, design.paths[index])
        rep_path = os.path.join(design.directory, design.paths[rep_index])
        is_rep = (index == rep_index)
        self.selected_model = None

        # Search for scripts that can perform some action using the clicked
        # model.  Such scripts must have the `*.sho' suffix and may be located
        # anywhere from the directory containing the models to any directory
        # below that.  Any scripts that are found will be used to populate a
        # drop-down menu.  If selected, the script will be called with sh as
        # the interpreter and the path to the model as the singular argument.

        directory = os.path.abspath(design.directory)
        sho_scripts = []

        while directory != os.path.abspath('/'):
            sho_pattern = os.path.join(directory, '*.sho')
            sho_scripts += glob.glob(sho_pattern)
            directory = os.path.dirname(directory)

        # Create and display the drop-down menu.

        file_menu = gtk.Menu()

        for script in sho_scripts:
            title = os.path.basename(os.path.splitext(script)[0])
            title = title[0].upper() + title[1:]
            title = title.replace('_', ' ')

            item = gtk.MenuItem(title)
            item.connect('activate',
                    lambda *args: try_to_run_command([script, path, rep_path]))
            file_menu.append(item)

        view_in_pymol = gtk.MenuItem("View model in pymol")
        view_in_pymol.connect('activate',
                lambda *args: try_to_run_command(['pymol', path]))
        file_menu.append(view_in_pymol)

        view_in_chimera = gtk.MenuItem("View model in chimera")
        view_in_chimera.connect('activate',
                lambda *args: try_to_run_command(['chimera', path]))
        file_menu.append(view_in_chimera)

        file_menu.append(gtk.SeparatorMenuItem())

        copy_path = gtk.MenuItem("Copy path to model")
        copy_path.connect('activate', self.on_copy_model_path, path)
        file_menu.append(copy_path)

        if index == design.representative:
            choose_rep = gtk.MenuItem("Reset representative")
            choose_rep.connect(
                'activate', self.on_set_representative, design, None)
        else:
            choose_rep = gtk.MenuItem("Set as representative")
            choose_rep.connect(
                'activate', self.on_set_representative, design, index)
        file_menu.append(choose_rep)

        file_menu.foreach(lambda item: item.show())
        ### Added a 'none'
        file_menu.popup(None, None, None, None, event.button, event.time)

    def on_copy_model_path(self, widget, path):
        import subprocess
        xsel = subprocess.Popen(['xsel', '-pi'], stdin=subprocess.PIPE)
        xsel.communicate(path)

    def on_set_representative(self, widget, design, index):
        design.representative = index
        self.update_plot()

    def on_edit_annotation(self, buffer):
        assert len(self.keys) == 1
        design = self.designs[self.keys[0]]
        bounds = buffer.get_bounds()
        design.notes = buffer.get_text(*bounds)

    def on_change_x_metric(self, widget):
        self.x_metric = widget.get_active_text()
        self.update_plot()

    def on_change_y_metric(self, widget):
        self.y_metric = widget.get_active_text()
        self.update_plot()
    
    def on_toggle_model_list(self, widget):
        if widget.get_active():
            self.show_model_list()
        else:
            self.hide_model_list()

    def on_toggle_filter_pane(self, widget):
        if widget.get_active():
            self.show_filter_pane()
        else:
            self.hide_filter_pane()

    def on_toggle_annotation_pane(self, widget):
        if widget.get_active():
            self.show_annotation_pane()
        else:
            self.hide_annotation_pane()

    def on_toggle_legend(self, widget):
        if widget.get_active():
            self.show_legend()
        else:
            self.hide_legend()

    def on_toggle_representative(self, widget):
        if widget.get_active():
            self.show_representative()
        else:
            self.hide_representative()

    def on_toggle_model_count(self, widget):
        if widget.get_active():
            self.show_model_count()
        else:
            self.hide_model_count()


    def normal_mode(self):
        self.set_focus(None)

        if self.toolbar._active == 'PAN':
            self.toolbar.pan()

        if self.toolbar._active == 'ZOOM':
            self.toolbar.zoom()

    def insert_mode(self):
        self.set_focus(self.notes)

    def search_mode(self):
        self.set_focus(self.search_form)

    def zoom_mode(self):
        self.toolbar.zoom()

    def pan_mode(self):
        self.toolbar.pan()

    def refocus_plot(self):
        self.toolbar.home()
        self.normal_mode()

    def next_design(self):
        selection = self.view.get_selection()
        model, paths = selection.get_selected_rows()
        num_paths = model.iter_n_children(None)
        if paths[-1][0] < model.iter_n_children(None) - 1:
            for path in paths: selection.unselect_path(path)
            selection.select_path(paths[-1][0] + 1)
            self.view.scroll_to_cell(paths[-1][0] + 1)

    def previous_design(self):
        selection = self.view.get_selection()
        model, paths = selection.get_selected_rows()
        if paths[0][0] > 0:
            for path in paths: selection.unselect_path(path)
            selection.select_path(paths[0][0] - 1)
            self.view.scroll_to_cell(paths[0][0] - 1)

    def cycle_x_metric(self, step=1):
        i = self.sorted_metrics.index(self.x_metric)
        i = (i + step) % len(self.sorted_metrics)
        if self.sorted_metrics[i] == self.y_metric:
            i = (i + step) % len(self.sorted_metrics)

        # Change the axis by programmatically selecting a new entry in the
        # corresponding drop-down menu in the toolbar.  This is incredibly
        # roundabout (and it kinda breaks encapsulation, although I consider
        # ModelViewer and ModelToolbar to be friends), but it's the only way I
        # know to keep the drop-down menu in sync.

        self.toolbar.x_axis_menu.set_active(i)

    def cycle_y_metric(self, step=1):
        i = self.sorted_metrics.index(self.y_metric)
        i = (i + step) % len(self.sorted_metrics)
        if self.sorted_metrics[i] == self.x_metric:
            i = (i + step) % len(self.sorted_metrics)

        # Change the axis by programmatically selecting a new entry in the
        # corresponding drop-down menu in the toolbar.  This is incredibly
        # roundabout (and it kinda breaks encapsulation, although I consider
        # ModelViewer and ModelToolbar to be friends), but it's the only way I
        # know to keep the drop-down menu in sync.

        self.toolbar.y_axis_menu.set_active(i)

    def reverse_cycle_x_metric(self):
        self.cycle_x_metric(-1)

    def reverse_cycle_y_metric(self):
        self.cycle_y_metric(-1)

    def save_selected_paths(self):
        chooser = gtk.FileChooserDialog(
                "Save selected paths",
                parent=self,
                action=gtk.FILE_CHOOSER_ACTION_SAVE,
                buttons=(
                    gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL,
                    gtk.STOCK_SAVE, gtk.RESPONSE_OK))

        chooser.set_current_folder(os.getcwd())
        chooser.set_current_name('selected_paths.txt')

        response = chooser.run()

        if response == gtk.RESPONSE_OK:
            selected_designs = [self.designs[key] for key in self.keys]
            with open(chooser.get_filename(), 'w') as file:
                file.writelines(
                        os.path.join(
                            design.directory,
                            design.paths[design.representative]) + '\n'
                        for design in selected_designs)

        chooser.destroy()

    def save_selected_funnels(self):
        from matplotlib.backends.backend_pdf import PdfPages
        import matplotlib.pyplot as plt

        selected_designs = [self.designs[key] for key in self.keys]

        chooser = gtk.FileChooserDialog(
                "Save selected funnels",
                parent=self,
                action=gtk.FILE_CHOOSER_ACTION_SAVE,
                buttons=(
                    gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL,
                    gtk.STOCK_SAVE, gtk.RESPONSE_OK))

        chooser.set_current_folder(os.getcwd())
        chooser.set_current_name('selected_funnels.pdf')

        response = chooser.run()

        if response == gtk.RESPONSE_OK:
            pdf = PdfPages(chooser.get_filename())

            for index, design in enumerate(selected_designs):
                plt.figure(figsize=(8.5, 11))
                plt.suptitle(design.directory)

                self.plot_models(plt.gca(), [design])

                pdf.savefig()
                plt.close()

            pdf.close()

        chooser.destroy()

    def hide_model_list(self):
        self.model_list.hide()
        self.model_list_toggle.set_active(False)
        if not self.filter_pane.props.visible:
            self.sidebar.hide()

    def show_model_list(self):
        self.model_list.show()
        self.model_list_toggle.set_active(True)
        self.sidebar.show()

    def toggle_model_list(self):
        if self.model_list.props.visible:
            self.hide_model_list()
        else:
            self.show_model_list()

    def hide_filter_pane(self):
        self.filter_pane.hide()
        self.filter_pane_toggle.set_active(False)
        if not self.model_list.props.visible:
            self.sidebar.hide()

    def show_filter_pane(self):
        self.filter_pane.show()
        self.filter_pane_toggle.set_active(True)
        self.sidebar.show()

    def toggle_filter_pane(self):
        if self.filter_pane.props.visible:
            self.hide_filter_pane()
        else:
            self.show_filter_pane()

    def hide_annotation_pane(self):
        self.annotation_pane.hide()
        self.annotation_pane_toggle.set_active(False)

    def show_annotation_pane(self):
        self.annotation_pane.show()
        self.annotation_pane_toggle.set_active(True)

    def toggle_annotation_pane(self):
        if self.annotation_pane.props.visible:
            self.hide_annotation_pane()
        else:
            self.show_annotation_pane()

    def hide_legend(self):
        if self.is_legend_visible:
            self.is_legend_visible = False
            self.legend_toggle.set_active(False)
            self.update_plot()

    def show_legend(self):
        if not self.is_legend_visible:
            self.is_legend_visible = True
            self.legend_toggle.set_active(True)
            self.update_plot()

    def hide_representative(self):
        if self.is_representative_visible:
            self.is_representative_visible = False
            self.representative_toggle.set_active(False)
            self.update_plot()

    def show_representative(self):
        if not self.is_representative_visible:
            self.is_representative_visible = True
            self.representative_toggle.set_active(True)
            self.update_plot()

    def toggle_legend(self):
        if self.is_legend_visible:
            self.hide_legend()
        else:
            self.show_legend()

    def hide_model_count(self):
        if self.is_model_count_visible:
            self.is_model_count_visible = False
            self.model_count_toggle.set_active(False)
            self.update_plot()

    def show_model_count(self):
        if not self.is_model_count_visible:
            self.is_model_count_visible = True
            self.model_count_toggle.set_active(True)
            self.update_plot()

    def toggle_model_count(self):
        if self.is_model_count_visible:
            self.hide_model_count()
        else:
            self.show_model_count()

    def plot_models(self, axes, designs, **kwargs):
        from itertools import count

        labels = kwargs.get('labels', None)
        x_metric = kwargs.get('x_metric', self.x_metric)
        y_metric = kwargs.get('y_metric', self.y_metric)

        # Define the colors that the plot will use.

        red =    '#ef2929', '#cc0000', '#a40000'
        orange = '#fcaf3e', '#f57900', '#ce5c00'
        yellow = '#fce94f', '#edd400', '#c4a000'
        green =  '#8ae234', '#73d216', '#4e9a06'
        blue =   '#729fcf', '#3465a4', '#204a87'
        purple = '#ad7fa8', '#75507b', '#5c3566'
        brown =  '#e9b96e', '#c17d11', '#8f5902'
        grey =   '#2e3436', '#555753', '#888a85', '#babdb6', '#d3d7cf', '#eeeeec'

        def color_from_cycle(index): #
            cycle = (blue[1], red[1], green[2], orange[1], purple[1], brown[1],
                     blue[0], red[0], green[1], orange[0], purple[0], brown[0])
            return cycle[index % len(cycle)]

        # Clear the axes and reset the axis labels

        axes.clear()
        axes.set_xlabel(self.metrics[x_metric].title)
        axes.set_ylabel(self.metrics[y_metric].title)

        # Plot the two axes.

        for index, design in enumerate(designs):
            rep = design.representative
            color = color_from_cycle(index)
            label = labels[index] if labels is not None else ''
            action = self.filter_pane.get_action()
            keep, drop = self.filter_pane.get_masks(design)

            x = design.get_metric(x_metric)
            y = design.get_metric(y_metric)

            # Scale the size of the points by the number of points.
            size = np.clip(7500 / max(len(x), 1), 2, 15)

            # Highlight the representative model.
            if self.is_representative_visible and keep[rep]:
                axes.scatter(
                        [x[rep]], [y[rep]],
                        s=60, c=yellow[1], marker='o', edgecolor='none',
                        label='_nolabel_')

            # Highlight the filtered points, if that's what the user wants.
            if action == 'highlight':
                axes.scatter(
                        x[drop], y[drop],
                        s=size, c=grey[4], marker='o', edgecolor='none',
                        label='_nolabel_')

            # Draw the whole score vs distance plot.
            lines = axes.scatter(
                    x[keep], y[keep],
                    s=size, c=color, marker='o', edgecolor='none',
                    label=label, picker=True)

            lines.paths = design.paths
            lines.design = design

        # Pick the axis limits based on the range of every design.  This is done
        # so you can scroll though every design without the axes changing size.

        def get_metric_limits(metric): #
            values = np.concatenate([x.get_metric(metric) for x in self])
            return self.metrics[metric].limits(values)

        x_min, x_max = get_metric_limits(x_metric)
        y_min, y_max = get_metric_limits(y_metric)

        x_pad = 0.05 * (x_max - x_min)
        y_pad = 0.05 * (y_max - y_min)

        axes.set_ylim(
            bottom=y_min - y_pad,
            top=y_max + y_pad,
        )
        axes.set_xlim(
            left=x_min - x_pad,
            right=x_max + x_pad,
        )

        # Draw guides for axes the that have them.

        x_guide = self.metrics[self.x_metric].guide
        y_guide = self.metrics[self.y_metric].guide

        if x_guide is not None:
            axes.axvline(x_guide, color=grey[3], linestyle='--')
        if y_guide is not None:
            axes.axhline(y_guide, color=grey[3], linestyle='--')

        # Draw the legend if the user enabled it.

        if self.is_legend_visible:
            axes.legend(loc='upper right')

        if self.is_model_count_visible:
            axes.annotate(
                    ', '.join(str(len(x)) for x in designs),
                    xy=(0, 1), xycoords='axes fraction',
                    xytext=(8, -8), textcoords='offset points',
                    verticalalignment='top',
            )


    def update_everything(self):
        self.update_annotations()
        self.update_plot()
        self.update_designs()

    def update_plot(self):
        designs = [self.designs[k] for k in self.keys]
        self.plot_models(self.axes, designs, labels=self.keys)
        self.canvas.draw()

    def update_annotations(self):
        if len(self.keys) == 1:
            design = self.designs[self.keys[0]]
            self.notes.get_buffer().set_text(design.notes)
            self.notes.set_sensitive(True)
        else:
            self.notes.set_sensitive(False)

    def update_designs(self):
        model = self.view.get_model()
        selector = self.view.get_selection()
        model.clear()

        def query_matches_design(design):
            needle = self.search_form.get_text()
            haystack = design.notes

            if needle.islower():
                haystack = haystack.lower()

            return needle in haystack

        for key in sorted(self.designs):
            if query_matches_design(self.designs[key]):
                model.append([key])

        selector.select_path((0,))


class ShowMyViolins (gtk.Window):

    def __init__(self, designs):

        # Setup the parent class.

        gtk.Window.__init__(self)
        self.add_events(gdk.EventMask.KEY_PRESS_MASK)
        self.connect('key-press-event', self.on_hotkey_press)
        self.set_icon_from_file(os.path.join(
            os.path.dirname(os.path.realpath(__file__)), 'icon.png'))

        # Setup the data members.

        self.designs = designs

        self.keys = list()
        self.selected_model = None
        self.is_legend_visible = False
        self.is_representative_visible = False
        self.is_model_count_visible = False

        self.metrics = {
                k: next(iter(self)).metrics[k] \
                for k in set.intersection(*[set(x.metrics) for x in self])
        }
        print(self.metrics)
        self.sorted_metrics = sorted(
                self.metrics,
                key=lambda k: (self.metrics[k].order, self.metrics[k].title)
        )
        print(self.sorted_metrics)
        self.x_metric = (
                default_x_metric
                if default_x_metric in self.metrics
                else self.sorted_metrics[0])
        # self.y_metric = (
        #         default_y_metric
        #         if default_y_metric in self.metrics
        #         else self.sorted_metrics[1])

        # Setup the GUI.

        self.connect('destroy', lambda x: gtk.main_quit())
        self.set_default_size(int(1.618 * 630), 630)

        menu_bar = self.setup_menu_bar()
        model_viewer = self.setup_model_viewer()

        vbox = gtk.VBox()
        vbox.pack_start(menu_bar, False)
        vbox.pack_end(model_viewer, True)

        self.add(vbox)
        self.set_border_width(3)
        self.update_everything()
        self.show_all()
        self.set_focus(None)

        n = len(self.designs)

        self.hide_model_list() if n == 1 else self.show_model_list()
        self.hide_filter_pane()
        self.show_annotation_pane()
        self.hide_legend()
        self.show_representative()
        self.hide_model_count()

    def __iter__(self):
        return iter(list(self.designs.values()))


    def setup_menu_bar(self):
        bar = gtk.MenuBar()

        # The "File" menu:
        menu = gtk.Menu()
        item = gtk.MenuItem("_File")
        item.set_submenu(menu)
        bar.append(item)

        item = gtk.MenuItem("Save selected paths")
        item.connect('activate', lambda _: self.save_selected_paths())
        menu.append(item)

        item = gtk.MenuItem("Save selected funnels")
        item.connect('activate', lambda _: self.save_selected_funnels())
        menu.append(item)

        # The "View" menu:
        menu = gtk.Menu()
        item = gtk.MenuItem("_View")
        item.set_submenu(menu)
        bar.append(item)

        item = self.model_list_toggle = gtk.CheckMenuItem("Model list")
        item.connect('activate', self.on_toggle_model_list)
        menu.append(item)

        item = self.filter_pane_toggle = gtk.CheckMenuItem("Filter pane")
        item.connect('activate', self.on_toggle_filter_pane)
        menu.append(item)

        item = self.annotation_pane_toggle = gtk.CheckMenuItem("Annotation pane")
        item.connect('activate', self.on_toggle_annotation_pane)
        menu.append(item)

        item = gtk.SeparatorMenuItem()
        menu.append(item)

        item = self.legend_toggle = gtk.CheckMenuItem("Legend")
        item.connect('activate', self.on_toggle_legend)
        menu.append(item)

        item = self.representative_toggle = gtk.CheckMenuItem("Representative")
        item.connect('activate', self.on_toggle_representative)
        menu.append(item)

        item = self.model_count_toggle = gtk.CheckMenuItem("Model count")
        item.connect('activate', self.on_toggle_model_count)
        menu.append(item)

        return bar

    def setup_model_viewer(self):
        plot = self.setup_plot()
        self.model_list = self.setup_model_list()
        self.filter_pane = self.setup_filter_pane()
        self.annotation_pane = self.setup_annotation_pane()

        self.sidebar = gtk.VPaned()
        self.sidebar.add1(self.model_list)
        self.sidebar.add2(self.filter_pane)

        bottombar = gtk.VPaned()
        bottombar.add1(plot)
        bottombar.add2(self.annotation_pane)

        viewer = gtk.HPaned()
        viewer.add1(self.sidebar)
        viewer.add2(bottombar)

        return viewer

    def setup_model_list(self):
        list_store = gtk.ListStore(str)

        text = gtk.CellRendererText()
        icon = gtk.CellRendererPixbuf()

        self.view = gtk.TreeView(list_store)
        self.view.set_model(list_store)
        self.view.set_rubber_banding(True)
        self.view.set_enable_search(False)
        self.view.set_headers_visible(False)

        columns = [
                ('Name', 'directory'),
        ]

        for index, parameters in enumerate(columns):
            title, attr = parameters

            def cell_data_func(column, cell, model, iter, attr): #
                key = model.get_value(iter, 0)
                design = self.designs[key]
                text = getattr(design, attr)
                cell.set_property('text', text)

            def sort_func(model, iter_1, iter_2, attr): #
                key_1 = model.get_value(iter_1, 0)
                key_2 = model.get_value(iter_2, 0)
                design_1 = self.designs[key_1]
                design_2 = self.designs[key_2]
                value_1 = getattr(design_1, attr)
                value_2 = getattr(design_2, attr)
                return cmp(value_1, value_2)

            list_store.set_sort_func(index, sort_func, attr);

            column = gtk.TreeViewColumn(title, text)
            column.set_cell_data_func(text, cell_data_func, attr)
            column.set_sort_column_id(index)
            self.view.append_column(column)

        selector = self.view.get_selection()
        selector.connect("changed", self.on_select_designs)
        ### If selection is weird, double check this!
        ### gtk.SelectionMode(3) should correspond to GTK_SELECTION_MULTIPLE
        selector.set_mode(gtk.SelectionMode(3))

        scroller = gtk.ScrolledWindow()
        scroller.set_policy(gtk.PolicyType.AUTOMATIC,
                gtk.PolicyType.AUTOMATIC)
        scroller.add(self.view)

        frame = gtk.Frame()
        frame.add(scroller)

        self.search_form = gtk.Entry()
        self.search_form.set_icon_from_stock(gtk.EntryIconPosition.SECONDARY, gtk.STOCK_FIND)

        search_buffer = self.search_form.get_buffer()
        search_buffer.connect('deleted-text', self.on_search_in_notes)
        search_buffer.connect('inserted-text', self.on_search_in_notes)

        vbox = gtk.VBox()
        vbox.pack_start(self.search_form, False, True, 0)
        vbox.pack_start(frame, True, True, 0)

        return vbox

    def setup_plot(self):
        figure = Figure(facecolor='#edecea')

        # Create the axes.

        self.axes = figure.add_axes((0.15, 0.15, 0.75, 0.75))
        self.axes.set_ylabel('Score')

        # Create the canvas.

        self.canvas = FigureCanvas(figure)
        self.canvas.mpl_connect('pick_event', self.on_select_model)
        self.canvas.mpl_connect('button_press_event', self.on_click_plot_mpl)
        self.canvas.mpl_connect('motion_notify_event', self.on_move_mouse_mpl)
        self.canvas.connect('button-press-event', self.on_click_plot_gtk)
        self.canvas.set_size_request(-1, 350)

        # Create the tool bar.

        self.toolbar = NavigationToolbar(self.canvas, self)

        # Place all the widgets.

        self.mouse_position = gtk.Label("")

        table = gtk.Table(3, 5)
        table.attach(self.toolbar, 0, 1, 0, 3)
        table.attach(self.mouse_position, 3, 4, 1, 2, xoptions=0, yoptions=0, xpadding=3)

        vbox = gtk.VBox()
        vbox.pack_start(self.canvas, True, True, 0)
        vbox.pack_start(table, False, True, 0)

        return vbox

    def setup_filter_pane(self):
        pane = FilterPane(self)
        pane.connect('updated', lambda _: self.update_plot())
        return pane

    def setup_annotation_pane(self):
        self.notes = gtk.TextView()
        self.notes.set_wrap_mode(gtk.WRAP_WORD)
        self.notes.set_size_request(-1, 100)
        self.notes.set_left_margin(3)
        self.notes.set_right_margin(3)
        self.notes.set_pixels_above_lines(3)
        self.notes.set_pixels_below_lines(3)
        self.notes.set_cursor_visible(True)
        self.notes.get_buffer().connect('changed', self.on_edit_annotation)

        scroll_window = gtk.ScrolledWindow()
        scroll_window.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)
        scroll_window.add(self.notes)

        frame = gtk.Frame()
        frame.add(scroll_window)

        return frame

    def setup_metric_menu(self, callback=None, initial_choice=None):
        try: self.metric_store
        except AttributeError:
            self.metric_store = gtk.ListStore(str, str)

            for key in self.sorted_metrics:
                metric = self.metrics[key]
                self.metric_store.append([metric.name, metric.title])

        cell = gtk.CellRendererText()
        ### previousoy ComboBox
        menu = gtk.ComboBoxText(model=self.metric_store)
        menu.pack_start(cell, True)
        menu.add_attribute(cell, 'text', 1)

        menu.set_active(0)
        for i, metric in enumerate(self.sorted_metrics):
            if metric == initial_choice:
                menu.set_active(i)

        if callback:
            menu.connect('changed', callback)

        return menu


    def on_hotkey_press(self, widget, event):
        key = gdk.keyval_name(event.keyval).lower()
        if event.state & gdk.CONTROL_MASK: key = 'ctrl-' + key
        if event.state & gdk.SHIFT_MASK: key = 'shift-' + key

        hotkeys = {
                'escape': self.normal_mode,
        }
        normal_mode_hotkeys = {
                'i': self.insert_mode,     'a': self.insert_mode,
                'z': self.zoom_mode,
                'x': self.pan_mode,
                'c': self.refocus_plot,
                'tab': self.cycle_y_metric,
                'space': self.cycle_x_metric,
                'shift-tab': self.reverse_cycle_y_metric,
                'shift-space': self.reverse_cycle_x_metric,
        }
        multi_design_hotkeys = {
                'j': self.next_design,     'f': self.next_design,
                'k': self.previous_design, 'd': self.previous_design,
                'slash': self.search_mode,
        }

        keep_focus = (
                gtk.Entry,
                gtk.TextView,
                gtk.Button,
                gtk.ComboBox
        )
        if not isinstance(self.get_focus(), keep_focus):
            hotkeys.update(normal_mode_hotkeys)
            if len(self.designs) > 1:
                hotkeys.update(multi_design_hotkeys)

        if key in hotkeys:
            hotkeys[key]()
            return True

    def on_search_in_notes(self, entry_buffer, *_):
        self.update_designs()

    def on_select_designs(self, selection):
        new_keys = []
        old_keys = self.keys[:]
        self.keys = []
        model, paths = selection.get_selected_rows()

        for path in paths:
            iter = model.get_iter(path)
            key = model.get_value(iter, 0)
            new_keys.append(key)

        # Don't change the order of designs that were already selected.  The
        # order affects how the color of the design in the score vs rmsd plot,
        # and things get confusing if it changes.

        for key in old_keys:
            if key in new_keys:
                self.keys.append(key)

        for key in new_keys:
            if key not in self.keys:
                self.keys.append(key)

        # Rename the window based on the current selection.

        subtitle = ""
        if len(self.keys) == 1:
            subtitle = " ({})".format(self.keys[0])
        if len(self.keys) > 1:
            subtitle = " ({}, ...)".format(self.keys[0])

        self.set_title("Show My Designs" + subtitle)

        # This is an efficiency thing.  The 'J' and 'K' hotkeys works in two
        # steps: first unselect everything and then select the next row in
        # order.  Redrawing the plot is expensive, so it's worthwhile to skip
        # redrawing after that first step.

        if self.keys:
            self.update_plot()
            self.update_annotations()

    def on_select_model(self, event):
        self.selected_model = event.ind[0], event.artist.design

    def on_move_mouse_mpl(self, event):
        if event.xdata is None or event.ydata is None:
            # The data coordinates will be None only if the mouse is outside
            # the data area.
            self.mouse_position.set_text("")
        else:
            coord = '{:0.2f}, {:0.2f}'.format(event.xdata, event.ydata)
            self.mouse_position.set_text(coord)

    def on_click_plot_mpl(self, event):
        pass

    def on_click_plot_gtk(self, widget, event):
        # Ignore any event that isn't a right button click or a left button
        # click with the control key held down.

        is_right_click = \
                (event.button == 3) or \
                (event.button == 1 and event.get_state() & gdk.CONTROL_MASK)

        if not is_right_click: return
        if self.toolbar._active == 'PAN': return
        if self.toolbar._active == 'ZOOM': return
        if self.selected_model is None: return

        # Figure out which model was clicked.

        index, design = self.selected_model
        rep_index = design.representative
        path = os.path.join(design.directory, design.paths[index])
        rep_path = os.path.join(design.directory, design.paths[rep_index])
        is_rep = (index == rep_index)
        self.selected_model = None

        # Search for scripts that can perform some action using the clicked
        # model.  Such scripts must have the `*.sho' suffix and may be located
        # anywhere from the directory containing the models to any directory
        # below that.  Any scripts that are found will be used to populate a
        # drop-down menu.  If selected, the script will be called with sh as
        # the interpreter and the path to the model as the singular argument.

        directory = os.path.abspath(design.directory)
        sho_scripts = []

        while directory != os.path.abspath('/'):
            sho_pattern = os.path.join(directory, '*.sho')
            sho_scripts += glob.glob(sho_pattern)
            directory = os.path.dirname(directory)

        # Create and display the drop-down menu.

        file_menu = gtk.Menu()

        for script in sho_scripts:
            title = os.path.basename(os.path.splitext(script)[0])
            title = title[0].upper() + title[1:]
            title = title.replace('_', ' ')

            item = gtk.MenuItem(title)
            item.connect('activate',
                    lambda *args: try_to_run_command([script, path, rep_path]))
            file_menu.append(item)

        view_in_pymol = gtk.MenuItem("View model in pymol")
        view_in_pymol.connect('activate',
                lambda *args: try_to_run_command(['pymol', path]))
        file_menu.append(view_in_pymol)

        view_in_chimera = gtk.MenuItem("View model in chimera")
        view_in_chimera.connect('activate',
                lambda *args: try_to_run_command(['chimera', path]))
        file_menu.append(view_in_chimera)

        file_menu.append(gtk.SeparatorMenuItem())

        copy_path = gtk.MenuItem("Copy path to model")
        copy_path.connect('activate', self.on_copy_model_path, path)
        file_menu.append(copy_path)

        if index == design.representative:
            choose_rep = gtk.MenuItem("Reset representative")
            choose_rep.connect(
                'activate', self.on_set_representative, design, None)
        else:
            choose_rep = gtk.MenuItem("Set as representative")
            choose_rep.connect(
                'activate', self.on_set_representative, design, index)
        file_menu.append(choose_rep)

        file_menu.foreach(lambda item: item.show())
        ### Added a 'none'
        file_menu.popup(None, None, None, None, event.button, event.time)

    def on_copy_model_path(self, widget, path):
        import subprocess
        xsel = subprocess.Popen(['xsel', '-pi'], stdin=subprocess.PIPE)
        xsel.communicate(path)

    def on_set_representative(self, widget, design, index):
        design.representative = index
        self.update_plot()

    def on_edit_annotation(self, buffer):
        assert len(self.keys) == 1
        design = self.designs[self.keys[0]]
        bounds = buffer.get_bounds()
        design.notes = buffer.get_text(*bounds)

    def on_change_x_metric(self, widget):
        self.x_metric = widget.get_active_text()
        self.update_plot()

    def on_change_y_metric(self, widget):
        self.y_metric = widget.get_active_text()
        self.update_plot()
    
    def on_toggle_model_list(self, widget):
        if widget.get_active():
            self.show_model_list()
        else:
            self.hide_model_list()

    def on_toggle_filter_pane(self, widget):
        if widget.get_active():
            self.show_filter_pane()
        else:
            self.hide_filter_pane()

    def on_toggle_annotation_pane(self, widget):
        if widget.get_active():
            self.show_annotation_pane()
        else:
            self.hide_annotation_pane()

    def on_toggle_legend(self, widget):
        if widget.get_active():
            self.show_legend()
        else:
            self.hide_legend()

    def on_toggle_representative(self, widget):
        if widget.get_active():
            self.show_representative()
        else:
            self.hide_representative()

    def on_toggle_model_count(self, widget):
        if widget.get_active():
            self.show_model_count()
        else:
            self.hide_model_count()


    def normal_mode(self):
        self.set_focus(None)

        if self.toolbar._active == 'PAN':
            self.toolbar.pan()

        if self.toolbar._active == 'ZOOM':
            self.toolbar.zoom()

    def insert_mode(self):
        self.set_focus(self.notes)

    def search_mode(self):
        self.set_focus(self.search_form)

    def zoom_mode(self):
        self.toolbar.zoom()

    def pan_mode(self):
        self.toolbar.pan()

    def refocus_plot(self):
        self.toolbar.home()
        self.normal_mode()

    def next_design(self):
        selection = self.view.get_selection()
        model, paths = selection.get_selected_rows()
        num_paths = model.iter_n_children(None)
        if paths[-1][0] < model.iter_n_children(None) - 1:
            for path in paths: selection.unselect_path(path)
            selection.select_path(paths[-1][0] + 1)
            self.view.scroll_to_cell(paths[-1][0] + 1)

    def previous_design(self):
        selection = self.view.get_selection()
        model, paths = selection.get_selected_rows()
        if paths[0][0] > 0:
            for path in paths: selection.unselect_path(path)
            selection.select_path(paths[0][0] - 1)
            self.view.scroll_to_cell(paths[0][0] - 1)

    def cycle_x_metric(self, step=1):
        i = self.sorted_metrics.index(self.x_metric)
        i = (i + step) % len(self.sorted_metrics)
        if self.sorted_metrics[i] == self.y_metric:
            i = (i + step) % len(self.sorted_metrics)

        # Change the axis by programmatically selecting a new entry in the
        # corresponding drop-down menu in the toolbar.  This is incredibly
        # roundabout (and it kinda breaks encapsulation, although I consider
        # ModelViewer and ModelToolbar to be friends), but it's the only way I
        # know to keep the drop-down menu in sync.

        self.toolbar.x_axis_menu.set_active(i)

    def cycle_y_metric(self, step=1):
        i = self.sorted_metrics.index(self.y_metric)
        i = (i + step) % len(self.sorted_metrics)
        if self.sorted_metrics[i] == self.x_metric:
            i = (i + step) % len(self.sorted_metrics)

        # Change the axis by programmatically selecting a new entry in the
        # corresponding drop-down menu in the toolbar.  This is incredibly
        # roundabout (and it kinda breaks encapsulation, although I consider
        # ModelViewer and ModelToolbar to be friends), but it's the only way I
        # know to keep the drop-down menu in sync.

        self.toolbar.y_axis_menu.set_active(i)

    def reverse_cycle_x_metric(self):
        self.cycle_x_metric(-1)

    def reverse_cycle_y_metric(self):
        self.cycle_y_metric(-1)

    def save_selected_paths(self):
        chooser = gtk.FileChooserDialog(
                "Save selected paths",
                parent=self,
                action=gtk.FILE_CHOOSER_ACTION_SAVE,
                buttons=(
                    gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL,
                    gtk.STOCK_SAVE, gtk.RESPONSE_OK))

        chooser.set_current_folder(os.getcwd())
        chooser.set_current_name('selected_paths.txt')

        response = chooser.run()

        if response == gtk.RESPONSE_OK:
            selected_designs = [self.designs[key] for key in self.keys]
            with open(chooser.get_filename(), 'w') as file:
                file.writelines(
                        os.path.join(
                            design.directory,
                            design.paths[design.representative]) + '\n'
                        for design in selected_designs)

        chooser.destroy()

    def save_selected_funnels(self):
        from matplotlib.backends.backend_pdf import PdfPages
        import matplotlib.pyplot as plt

        selected_designs = [self.designs[key] for key in self.keys]

        chooser = gtk.FileChooserDialog(
                "Save selected funnels",
                parent=self,
                action=gtk.FILE_CHOOSER_ACTION_SAVE,
                buttons=(
                    gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL,
                    gtk.STOCK_SAVE, gtk.RESPONSE_OK))

        chooser.set_current_folder(os.getcwd())
        chooser.set_current_name('selected_funnels.pdf')

        response = chooser.run()

        if response == gtk.RESPONSE_OK:
            pdf = PdfPages(chooser.get_filename())

            for index, design in enumerate(selected_designs):
                plt.figure(figsize=(8.5, 11))
                plt.suptitle(design.directory)

                self.plot_models(plt.gca(), [design])

                pdf.savefig()
                plt.close()

            pdf.close()

        chooser.destroy()

    def hide_model_list(self):
        self.model_list.hide()
        self.model_list_toggle.set_active(False)
        if not self.filter_pane.props.visible:
            self.sidebar.hide()

    def show_model_list(self):
        self.model_list.show()
        self.model_list_toggle.set_active(True)
        self.sidebar.show()

    def toggle_model_list(self):
        if self.model_list.props.visible:
            self.hide_model_list()
        else:
            self.show_model_list()

    def hide_filter_pane(self):
        self.filter_pane.hide()
        self.filter_pane_toggle.set_active(False)
        if not self.model_list.props.visible:
            self.sidebar.hide()

    def show_filter_pane(self):
        self.filter_pane.show()
        self.filter_pane_toggle.set_active(True)
        self.sidebar.show()

    def toggle_filter_pane(self):
        if self.filter_pane.props.visible:
            self.hide_filter_pane()
        else:
            self.show_filter_pane()

    def hide_annotation_pane(self):
        self.annotation_pane.hide()
        self.annotation_pane_toggle.set_active(False)

    def show_annotation_pane(self):
        self.annotation_pane.show()
        self.annotation_pane_toggle.set_active(True)

    def toggle_annotation_pane(self):
        if self.annotation_pane.props.visible:
            self.hide_annotation_pane()
        else:
            self.show_annotation_pane()

    def hide_legend(self):
        if self.is_legend_visible:
            self.is_legend_visible = False
            self.legend_toggle.set_active(False)
            self.update_plot()

    def show_legend(self):
        if not self.is_legend_visible:
            self.is_legend_visible = True
            self.legend_toggle.set_active(True)
            self.update_plot()

    def hide_representative(self):
        if self.is_representative_visible:
            self.is_representative_visible = False
            self.representative_toggle.set_active(False)
            self.update_plot()

    def show_representative(self):
        if not self.is_representative_visible:
            self.is_representative_visible = True
            self.representative_toggle.set_active(True)
            self.update_plot()

    def toggle_legend(self):
        if self.is_legend_visible:
            self.hide_legend()
        else:
            self.show_legend()

    def hide_model_count(self):
        if self.is_model_count_visible:
            self.is_model_count_visible = False
            self.model_count_toggle.set_active(False)
            self.update_plot()

    def show_model_count(self):
        if not self.is_model_count_visible:
            self.is_model_count_visible = True
            self.model_count_toggle.set_active(True)
            self.update_plot()

    def toggle_model_count(self):
        if self.is_model_count_visible:
            self.hide_model_count()
        else:
            self.show_model_count()

    def plot_models(self, axes, designs, **kwargs):
        from itertools import count

        labels = kwargs.get('labels', None)
        x_metric = kwargs.get('x_metric', self.x_metric)
        # y_metric = kwargs.get('y_metric', self.y_metric)

        # Define the colors that the plot will use.

        red =    '#ef2929', '#cc0000', '#a40000'
        orange = '#fcaf3e', '#f57900', '#ce5c00'
        yellow = '#fce94f', '#edd400', '#c4a000'
        green =  '#8ae234', '#73d216', '#4e9a06'
        blue =   '#729fcf', '#3465a4', '#204a87'
        purple = '#ad7fa8', '#75507b', '#5c3566'
        brown =  '#e9b96e', '#c17d11', '#8f5902'
        grey =   '#2e3436', '#555753', '#888a85', '#babdb6', '#d3d7cf', '#eeeeec'

        def color_from_cycle(index): #
            cycle = (blue[1], red[1], green[2], orange[1], purple[1], brown[1],
                     blue[0], red[0], green[1], orange[0], purple[0], brown[0])
            return cycle[index % len(cycle)]

        # Clear the axes and reset the axis labels

        axes.clear()
        axes.set_xlabel(self.metrics[x_metric].title)
        

        all_designs = {}
        
        if designs:
            for index, design in enumerate(designs):
                this_design = design.get_metric(x_metric).to_frame()
                
                this_design['name'] = design.directory
                all_designs[design.directory] = this_design
            data = pd.concat(all_designs)

            lines = sns.violinplot(x=data[x_metric], y = data['name'], ax = axes, inner = 'point')

            

        # Pick the axis limits based on the range of every design.  This is done
        # so you can scroll though every design without the axes changing size.

        def get_metric_limits(metric): #
            values = np.concatenate([x.get_metric(metric) for x in self])
            return self.metrics[metric].limits(values)

        x_min, x_max = get_metric_limits(x_metric)
        # y_min, y_max = get_metric_limits(y_metric)

        x_pad = 0.05 * (x_max - x_min)
        # y_pad = 0.05 * (y_max - y_min)

        # axes.set_ylim(
        #     bottom=y_min - y_pad,
        #     top=y_max + y_pad,
        # )
        axes.set_xlim(
            left=x_min - x_pad,
            right=x_max + x_pad,
        )

        # Draw guides for axes the that have them.

        x_guide = self.metrics[self.x_metric].guide
        # y_guide = self.metrics[self.y_metric].guide

        if x_guide is not None:
            axes.axvline(x_guide, color=grey[3], linestyle='--')
        # if y_guide is not None:
        #     axes.axhline(y_guide, color=grey[3], linestyle='--')

        # Draw the legend if the user enabled it.

        if self.is_legend_visible:
            axes.legend(loc='upper right')

        if self.is_model_count_visible:
            axes.annotate(
                    ', '.join(str(len(x)) for x in designs),
                    xy=(0, 1), xycoords='axes fraction',
                    xytext=(8, -8), textcoords='offset points',
                    verticalalignment='top',
            )


    def update_everything(self):
        self.update_annotations()
        self.update_plot()
        self.update_designs()

    def update_plot(self):
        designs = [self.designs[k] for k in self.keys]
        self.plot_models(self.axes, designs, labels=self.keys)
        self.canvas.draw()

    def update_annotations(self):
        if len(self.keys) == 1:
            design = self.designs[self.keys[0]]
            self.notes.get_buffer().set_text(design.notes)
            self.notes.set_sensitive(True)
        else:
            self.notes.set_sensitive(False)

    def update_designs(self):
        model = self.view.get_model()
        selector = self.view.get_selection()
        model.clear()

        def query_matches_design(design):
            needle = self.search_form.get_text()
            haystack = design.notes

            if needle.islower():
                haystack = haystack.lower()

            return needle in haystack

        for key in sorted(self.designs):
            if query_matches_design(self.designs[key]):
                model.append([key])

        selector.select_path((0,))


class FigureCanvas (FigureCanvasGTK3Agg):

    def __init__(self, figure):
        FigureCanvasGTK3Agg.__init__(self, figure)

    def button_press_event(self, widget, event):
        FigureCanvasGTK3Agg.button_press_event(self, widget, event)
        return False


class NavigationToolbar (NavigationToolbar2GTK3):

    toolitems = (
        ('Home', 'Reset original view', 'home', 'home'),
        ('Back', 'Back to previous view', 'back', 'back'),
        ('Forward', 'Forward to next view', 'forward', 'forward'),
        #(None, None, None, None),
        ('Pan', 'Pan axes with left mouse, zoom with right', 'move', 'pan'),
        ('Zoom', 'Zoom to rectangle', 'zoom_to_rect', 'zoom'),
        #(None, None, None, None),
        ('Save', 'Save the figure', 'filesave', 'save_figure'),
    )

    def __init__(self, canvas, parent):
        NavigationToolbar2GTK3.__init__(self, canvas, parent)

        self.x_axis_menu = parent.setup_metric_menu(
                parent.on_change_x_metric, parent.x_metric)
        # self.y_axis_menu = parent.setup_metric_menu(
                # parent.on_change_y_metric, parent.y_metric)

        table = gtk.Table(3, 4)
        table.attach(gtk.SeparatorToolItem(), 0, 1, 0, 3)
        # table.attach(self.y_axis_menu, 1, 2, 1, 2, xoptions=0, yoptions=0)
        table.attach(gtk.Label(' vs. '), 2, 3, 1, 2, xoptions=0, yoptions=0)
        table.attach(self.x_axis_menu, 3, 4, 1, 2, xoptions=0, yoptions=0)

        tool_item = gtk.ToolItem()
        tool_item.add(table)

        self.insert(tool_item, len(self.toolitems))

    def set_message(self, message):
        pass


class FilterPane(gtk.Table):

    __gsignals__ = {
            'updated' : (
                gobject.SIGNAL_RUN_LAST,
                gobject.TYPE_NONE,
                (),
            )
    }

    def __init__(self, master):
        gtk.Table.__init__(self)
        self.master = master
        self.filters = []
        self.action_menu = self.make_action_menu()
        self.add_button = self.make_add_button()
        self.update_num_rows()

    def get_action(self):
        ### Another fix to get_active_text
        #model = self.action_menu.get_model()
        #text = model.get_value(
        #        self.action_menu.get_active_iter(),
        #        self.action_menu.props.entry_text_column
        #        )
        return self.action_menu.get_active_text().lower()
        #return text.lower()

    def get_masks(self, design):
        keep = np.ones(len(design), dtype='bool')

        for filter in self.filters:
            name = filter.get_name()
            op = filter.get_operator()
            try: threshold = float(filter.get_threshold())
            except: continue

            metric = design.get_metric(name)

            if op == '>': result = metric > threshold
            if op == '<': result = metric < threshold
            if op == '>=': result = metric >= threshold
            if op == '<=': result = metric <= threshold
            if op == '==': result = metric == threshold
            if op == '!=': result = metric != threshold

            filter.update_counter(sum(result), len(result))
            keep &= result

        return keep, np.logical_not(keep)

    def make_add_button(self):
        button = make_stock_button(gtk.STOCK_ADD)
        button.connect('clicked', lambda _: self.add_filter())
        align = gtk.Alignment(0.0, 0.5)
        align.add(button)
        return align

    def make_action_menu(self):
        ### Check here if filtering is messed up
        #combo_box = gtk.combo_box_new_text()
        combo_box = gtk.ComboBoxText()
        combo_box.append_text("Highlight")
        combo_box.append_text("Hide")
        combo_box.set_active(0)
        combo_box.connect('changed', lambda _: self.emit('updated'))
        return combo_box

    def add_filter(self):
        filter = self.Filter(self)
        self.filters.append(filter)
        self.update_num_rows()

    def remove_filter(self, filter):
        self.filters.remove(filter)
        self.update_num_rows()

    def update_num_rows(self):
        # Remove everything.
        for child in self.get_children():
            self.remove(child)

        # Re-create the labels.
        action_label = gtk.Label("Action:")
        action_label.set_alignment(0.0, 0.5)

        filter_label = gtk.Label("Filters:")
        filter_label.set_alignment(0.0, 0.5)

        # Make the table the right size.
        rows = max(len(self.filters) + 2, 2)
        self.resize(rows, 6)

        # Re-attach everything to the table.
        fill = dict(xoptions=gtk.FILL, yoptions=gtk.FILL)

        self.attach(action_label,       0, 1, 0, 1, **fill)
        self.attach(self.action_menu,   1, 2, 0, 1, **fill)
        self.attach(filter_label,       0, 1, 1, 2, **fill)

        i = 0
        for filter in self.filters:
            filter.attach(i+1, **fill)
            i += 1

        self.attach(self.add_button,    1, 2, i+1, i+2, **fill)
        self.emit('updated')
        self.show_all()


    class Filter(object):

        def __init__(self, table):
            self.table = table
            self.filter_menu = table.master.setup_metric_menu()
            self.operator_menu = make_operator_menu()
            self.threshold_entry = gtk.Entry()
            self.threshold_entry.set_width_chars(5)
            self.delete_button = make_stock_button(gtk.STOCK_CANCEL)
            self.counter = gtk.Label()

            self.filter_menu.connect('changed', lambda _: table.emit('updated'))
            self.operator_menu.connect('changed', lambda _: table.emit('updated'))
            self.threshold_entry.connect('activate', lambda _: table.emit('updated'))
            self.delete_button.connect('clicked', lambda _: table.remove_filter(self))

        def __repr__(self):
            return '<Filter "{} {} {}">'.format(
                    self.filter_menu.get_active_text(),
                    self.operator_menu.get_active_text(),
                    self.threshold_entry.get_text() or '???')

        def get_name(self):
            return self.filter_menu.get_active_text()

        def get_operator(self):
            return self.operator_menu.get_active_text()

        def get_threshold(self):
            return self.threshold_entry.get_text()

        def attach(self, i, **fill):
            self.table.attach(self.filter_menu,      1, 2, i, i+1, **fill)
            self.table.attach(self.operator_menu,    2, 3, i, i+1, **fill)
            self.table.attach(self.threshold_entry,  3, 4, i, i+1, **fill)
            self.table.attach(self.delete_button,    4, 5, i, i+1, **fill)
            self.table.attach(self.counter,          5, 6, i, i+1, **fill)

        def update_counter(self, num_kept, num_total):
            self.counter.set_text('{}/{}'.format(num_kept, num_total))



class MetricInfo(object):

    def __init__(self, name, title, order, guide, limits):
        self.name = name
        self.title = title
        self.order = order
        self.guide = guide
        self.limits = limits

    def __repr__(self):
        return '<MetricInfo name="{0}">'.format(self.name)



def make_stock_button(stock):
    image = gtk.Image()
    image.set_from_stock(stock, gtk.IconSize.BUTTON)
    button = gtk.Button()
    button.add(image)
    return button

def make_operator_menu():
    ### ComboBoxText
    #combo_box = gtk.combo_box_new_text()
    combo_box = gtk.ComboBoxText()
    options = '<', '>', '<=', '>=', '=', '!='
    for option in options:
        combo_box.append_text(option)
    combo_box.set_active(0)
    return combo_box


default_x_metric = 'restraint_dist'
default_y_metric = 'total_score'

metric_titles = {
        'total_score': 'Total Score (REU)',
        'loop_rmsd': 'Loop RMSD (Å)',
        'delta_buried_unsats': 'Δ Buried Unsats',
}

metric_orders = {
}

metric_guides = {
        'loop_rmsd': 1.0,
}

metric_limits = {
        'total_score': lambda x: (
            min(x),
            np.percentile(x, 85)),

        'loop_rmsd': lambda x: (
            0.025 * max(x),
            max(x)),
}


def get_metric_title(metric, design=None):
    naive_title = metric.replace('_', ' ').title()
    return metric_titles.get(metric, naive_title)

def get_metric_order(metric, design=None):
    return metric_orders.get(metric)

def get_metric_guide(metric, design=None):
    return metric_guides.get(metric)

def get_metric_limits(metric, design=None):
    return metric_limits.get(metric, lambda x: (min(x), max(x)))


def show_my_designs(directories, use_cache=True, launch_gui=True, fork_gui=True):
    try:
        designs = load_designs(directories, use_cache=use_cache)
        #designs['relax_holo/01_relax_models/outputs/']._load_metrics

        if designs and launch_gui:
            # If the user wants to run in a background process, try to fork.
            # But for some reason fork() doesn't seem to work on Macs, so just
            # run the GUI in the main process if anything goes wrong.
            try:
                if fork_gui and os.fork():
                    sys.exit()
            except Exception:
                pass

            gui = ShowMyDesigns(designs)
            gtk.main()

    except KeyboardInterrupt:
        print()


def show_my_violins(directories, use_cache=True, launch_gui=True, fork_gui=True):
    try:
        designs = load_designs(directories, use_cache=use_cache)
        #designs['relax_holo/01_relax_models/outputs/']._load_metrics

        if designs and launch_gui:
            # If the user wants to run in a background process, try to fork.
            # But for some reason fork() doesn't seem to work on Macs, so just
            # run the GUI in the main process if anything goes wrong.
            try:
                if fork_gui and os.fork():
                    sys.exit()
            except Exception:
                pass

            gui = ShowMyViolins(designs)
            print("designs", designs)
            gtk.main()

    except KeyboardInterrupt:
        print()


def load_designs(directories, use_cache=True):
    designs = collections.OrderedDict()

    for directory in directories:
        try:
            designs[directory] = Design(directory, use_cache)
        except IOError as error:
            if str(error):
                print("Error:", str(error))
            else:
                raise

    return designs

def parse_records_from_pdbs(pdb_paths):
    records = []

    for i, path in enumerate(pdb_paths):

        # Update the user on our progress, because this is often slow.

        sys.stdout.write("\rReading '{}' [{}/{}]".format(
            os.path.dirname(path), i+1, len(pdb_paths)))
        sys.stdout.flush()

        # Read the PDB file, which we are assuming is gzipped.

        try:
            def smart_open(path): #
                if path.endswith('.gz'): return gzip.open(path)
                else: return open(path)

            with smart_open(path) as file:
                lines = file.readlines()

        except IOError:
            print("\nFailed to read '{}'".format(path))
            continue

        # Parse the pdb file.  This method may be overloaded to parse
        # different kinds of information.

        record = {'path': os.path.basename(path)}
        parse_record_from_pdb(record, path, lines)
        records.append(record)

    if pdb_paths: print()
    return records

def parse_record_from_pdb(record, pdb_path, lines):
    # Get different information from different lines in the PDB file.  Some
    # of these lines are specific to certain simulations.

    for line in lines:
        if line.startswith('total_score'):
            record['total_score'] = float(line.split()[1])

        if line.startswith('pose'):
            record['total_score'] = float(line.split()[-1])

        if line.startswith('loop_backbone_rmsd'):
            record['loop_rmsd'] = float(line.split()[1])

        if line.startswith('delta_buried_unsats'):
            record['delta_buried_unsats'] = float(line.split()[1])

def try_to_run_command(command):
    with open(os.devnull, 'w') as devnull:
        try: subprocess.Popen(command, stdout=devnull)
        except OSError as error:
            message = gtk.MessageDialog(
                    parent=None,
                    flags=0,
                    type=gtk.MESSAGE_ERROR,
                    buttons=gtk.BUTTONS_OK,
            )
            message.set_markup("<b>Failed to run {}</b>".format(command[0]))
            message.format_secondary_text(str(error))
            message.run()
            message.destroy()


def main():
    import docopt
    args = docopt.docopt(__doc__)

    if args['--version']:
        from show_my_designs import __version__
        print('{0} {1} (python {2[0]}.{2[1]})'.format(
                os.path.basename(sys.argv[0]), __version__, sys.version_info))
        raise SystemExit

    show_my_designs(
            args['<pdb_directories>'],
            use_cache=not args['--force'],
            launch_gui=not args['--quiet'],
            fork_gui=not args['--no-fork'],
    )
