"""
TrajectoryImageGenerator: A class for generating images and movies from molecular dynamics trajectories.
This class provides functionality to create images from molecular dynamics trajectories using NGLView,
and compile those images into a movie. It supports parallel image generation using multiple threads.

Copyright (C) 2025 Emile de Bruyn

This library is free software; you can redistribute it and/or modify it under the terms of the
GNU Lesser General Public License as published by the Free Software Foundation; either version 2.1
of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with this library;
if not, see <https://www.gnu.org/licenses/>.

Recent Changes:
---------------
- Added smoothing functionality for trajectory positions using a Gaussian-like filter.
- Improved image generation by ensuring proper handling of frame indices and rendering.
- Enhanced thread management with safe thread stopping and cleanup.
- Updated movie creation to use a predefined list of frame indices for consistent output.

Example usage:
--------------
# Initialize the image generator with a trajectory, output folder, indices, and number of views
image_generator = TrajectoryImageGenerator(
    smoothed_traj, 
    '7wii_frames_aa', 
    range(0, int(len(smoothed_traj.trajectory)/10), 1), 
    num_views=5
)

# Generate images in parallel
image_generator.run()

# Wait until the image generation is completed, then create a movie
image_generator.make_movie('7wii_frames_aa.mp4', fps=30)

# Cleanup resources and close views
image_generator.cleanup()

Attributes:
-----------
traj : MDAnalysis.Universe
    The molecular dynamics trajectory to visualize.
image_folder : str
    The folder where generated images will be saved.
indices : list
    The list of frame indices to generate images for.
smoothing_window : int
    The size of the window used for smoothing trajectory positions.
num_views : int
    The number of NGLView widgets (views) to use for parallel image generation.
views : list
    A list of NGLView widgets created for the trajectory.
progress_bar : tqdm.notebook.tqdm
    A progress bar to track the image generation process.
threads : list
    A list of threads used for parallel image generation.
rep_args : list
    A list of representation arguments for customizing NGLView visualizations.
selection : str
    A selection string for specifying which atoms to visualize.

Methods:
--------
get_smooth_positions():
    Smooth the trajectory positions using a Gaussian-like filter.
smooth_trajectory():
    Apply smoothing to the trajectory positions.
create_views():
    Create NGLView widgets for the trajectory.
base64_to_ndarray(value):
    Convert base64 image data to a numpy array.
get_image_from_view(view, frame_index):
    Render an image from the view at a specific frame.
generate_images(view, subset_indices):
    Generate images for a subset of indices.
run():
    Start parallel image generation.
_async_raise(tid, exctype):
    Raise an exception in a thread with a specific ID.
_stop_thread(thread):
    Stop a thread by raising SystemExit in it.
sorted_alphanumeric(list):
    Sort a list in alphanumeric order.
cleanup():
    Cleanup the progress bar and close all views.
make_movie(output_file, fps=30):
    Create a movie from the generated images.
"""

import os
import threading
import time
import base64
import io
import tqdm
import nglview as nv
import numpy as np
import re
from moviepy.video.io.ImageSequenceClip import ImageSequenceClip
import ctypes

import matplotlib.image as mpimg


class TrajectoryImageGenerator:
    def __init__(self, traj, image_folder, indices, smoothing_window=10, num_views=10, rep_args=None, selection=None):
        self.traj = traj
        self.smoothing_window = smoothing_window
        self.image_folder = image_folder
        self.indices = indices
        self.views = []
        self.progress_bar = None
        self.num_views = num_views
        self.rep_args = rep_args
        self.selection = selection
        self.threads = []
        
        # Smooth the trajectory
        if self.smoothing_window > 1:
            self.smooth_trajectory()

        # Ensure the image folder exists
        os.makedirs(self.image_folder, exist_ok=True)

        # Create views upon initialisation that are used in a new cell to generate images
        self.create_views()

    @staticmethod
    def get_smooth_positions(self):
        """Smooth the trajectory using a Gaussian filter."""
        # Implement smoothing logic here if needed
        smoothed_positions = np.zeros_like(self.traj.trajectory.timeseries())
        for i in range(len(self.traj.trajectory)):
            start = max(0, i - self.smoothing_window // 2)
            end = min(len(self.traj.trajectory), i + self.smoothing_window // 2 + 1)
            # smoothed_positions[i] = np.mean(traj.trajectory[start:end].positions, axis=0)
            smoothed_positions[:, i, :] = np.mean(self.traj.trajectory.timeseries()[:, start:end, :], axis=1)
        return smoothed_positions
    
    def smooth_trajectory(self):
        smoothed_positions = self.get_smooth_positions(self)

        for i, ts in enumerate(self.traj.trajectory):
            ts.positions = smoothed_positions[:, i, :]

    def create_views(self):
        """Create NGLView widgets for the trajectory."""
        for i in range(self.num_views):
            view = nv.show_mdanalysis(self.traj)
            view.clear_representations()
            if self.rep_args is not None:
                for rep_arg in self.rep_args:
                    view.add_representation(**rep_arg)
            if self.selection is None:
                view.add_ball_and_stick(selection='protein')
                view.add_ball_and_stick(selection='nucleic')
            else:
                view.add_ball_and_stick(selection=self.selection)
            view.layout = {'width': '100%', 'height': '600px'}
            view.handle_resize()
            view.center()
            view.control.zoom(0.3)
            self.views.append(view)
            display(view)

    @staticmethod
    def base64_to_ndarray(value):
        """Convert base64 image data to a numpy array."""
        im_bytes = base64.b64decode(value)
        im_bytes = io.BytesIO(im_bytes)
        return mpimg.imread(im_bytes, format='PNG')

    def get_image_from_view(self, view, frame_index):
        """Render an image from the view at a specific frame."""
        view.frame = frame_index
        img = view.render_image(frame=frame_index, factor=4, transparent=True)
        while not img.value:
            time.sleep(0.1)
        im_data = view._image_data
        return self.base64_to_ndarray(im_data)

    def generate_images(self, view, subset_indices):
        """Generate images for a subset of indices."""
        for i in subset_indices:
            img_path = os.path.join(self.image_folder, f'frame_{i:06d}.png')
            mpimg.imsave(img_path, self.get_image_from_view(view, i))
            self.progress_bar.update(1)

    def run(self):
        """Create generate images in parallel. Use in a new cell."""
        indices_split = np.array_split(self.indices, len(self.views))
        indices_split = [subset.tolist() for subset in indices_split]

        # Create a progress bar
        self.progress_bar = tqdm.notebook.tqdm(total=len(self.indices), desc="Generating images", position=0, leave=True)

        for view, subset_indices in zip(self.views, indices_split):
            thread = threading.Thread(
                target=self.generate_images,
                args=(view, subset_indices),
                name="generate_images",
                daemon=True
            )
            self.threads.append(thread)
            thread.start()
   
    @staticmethod
    def _async_raise(tid, exctype):
        """Raises an exception in the threads with id tid"""
        res = ctypes.pythonapi.PyThreadState_SetAsyncExc(ctypes.c_long(tid), ctypes.py_object(exctype))
        if res == 0:
            raise ValueError("Invalid thread id")
        elif res != 1:
            # If it returns a number greater than one, we're in trouble, and we call it again with 0 to revert.
            ctypes.pythonapi.PyThreadState_SetAsyncExc(ctypes.c_long(tid), None)
            raise SystemError("PyThreadState_SetAsyncExc failed")

    @staticmethod
    def _stop_thread(self, thread):
        """Stops a thread by raising SystemExit in it."""
        self._async_raise(thread.ident, SystemExit)

    @staticmethod
    def sorted_alphanumeric(list):
        """Sorts a list in alphanumeric order."""
        convert = lambda text: int(text) if text.isdigit() else text.lower()
        alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
        return sorted(list, key=alphanum_key)

    def cleanup(self):
        """Cleanup the progress bar and close all views."""
        if self.progress_bar:
            self.progress_bar.close()
        for view in self.views:
            view.close()
        for thread in self.threads:
            try:
                self._stop_thread(thread)
            except TypeError:
                continue
 
    def make_movie(self, output_file, fps=30):
        """Create a movie from the generated images."""

        # Collect all image files in the folder
        # image_files = self.sorted_alphanumeric(
        #     [os.path.join(self.image_folder, f) for f in os.listdir(self.image_folder) if f.endswith('.png')]
        # )

        image_files = [os.path.join(self.image_folder, f'frame_{i:06d}.png') for i in self.indices]

        # Create a video clip from the image sequence
        clip = ImageSequenceClip(image_files, fps=fps)

        # Write the video file
        clip.write_videofile(output_file, codec='libx264')


# Example usage:
# >>> image_generator = TrajectoryImageGenerator(smoothed_traj, '7wii_frames_aa', range(0, int(len(smoothed_traj.trajectory)/10), 1), num_views=5)

# New cell:
# >>> image_generator.run()

# Wait until completed, then in a new cell:
# >>> image_generator.make_movie('7wii_frames_aa.mp4', fps=30)
# >>> image_generator.cleanup()
