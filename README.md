# TrajectoryImageGenerator

`TrajectoryImageGenerator` is a Python class designed to generate high-quality images and movies from molecular dynamics (MD) trajectories. It leverages the power of [NGLView](https://github.com/nglviewer/nglview) for visualization and [MoviePy](https://zulko.github.io/moviepy/) for video creation. The class supports parallel image generation using multiple threads, making it efficient for large-scale trajectory visualization tasks.

---

## Features

- **Trajectory Visualization**: Create images of molecular dynamics trajectories using NGLView.
- **Movie Generation**: Compile generated images into a movie with customizable frame rates.
- **Smoothing**: Smooth trajectory positions using a Gaussian-like filter for better visualization.
- **Parallel Processing**: Generate images in parallel using multiple threads for faster processing.
- **Customizable Representations**: Add custom visual representations (e.g., ball-and-stick) for proteins, nucleic acids, or specific atom selections.
- **Progress Tracking**: Monitor the image generation process with a progress bar.

---

## Installation

Ensure you have the following dependencies installed:

- Python 3.7+
- [MDAnalysis](https://www.mdanalysis.org/)
- [NGLView](https://github.com/nglviewer/nglview)
- [MoviePy](https://zulko.github.io/moviepy/)
- [Matplotlib](https://matplotlib.org/)

Install the required Python packages using pip:

```bash
pip install mdanalysis nglview moviepy matplotlib tqdm
```

## Usage
### Example Workflow
 1. Initialize the Image Generator
    Create an instance of the `TrajectoryImageGenerator` class with your trajectory, output folder, and other parameters.
    ```Python
    from TrajectoryImageGenerator import TrajectoryImageGenerator

    image_generator = TrajectoryImageGenerator(
        traj=smoothed_traj, 
        image_folder='output_frames', 
        indices=range(0, len(smoothed_traj.trajectory), 10), 
        num_views=5
    )
    ```
 2. Generate Images
    ```Python
    image_generator.run()
    ```
 3. Create a Movie
    Once image generation is complete, compile the images into a movie.
    ```Python
    image_generator.make_movie('trajectory_movie.mp4', fps=30)
    ```
 4. Cleanup Resources
    Close all views and clean up resources.
    ```Python
    image_generator.cleanup()
    ```

## Class Details
### Attributes
 - `traj`: The MD trajectory to visualize (an MDAnalysis.Universe object).
 - `image_folder`: Directory where generated images will be saved.
 - `indices`: List of frame indices to generate images for.
 - `smoothing_window`: Size of the window used for smoothing trajectory positions.
 - `num_views`: Number of NGLView widgets (views) to use for parallel image generation.
 - `rep_args`: List of representation arguments for customizing visualizations.
 - `selection`: Atom selection string for visualization.

### Methods
 - `smooth_trajectory()`: Smooth trajectory positions using a Gaussian-like filter.
 - `create_views()`: Create NGLView widgets for trajectory visualization.
 - `generate_images(view, subset_indices)`: Generate images for a subset of frame indices.
 - `run()`: Start parallel image generation.
 - `make_movie(output_file, fps=30)`: Create a movie from the generated images.
 - `cleanup()`: Cleanup resources and close all views.

## Advanced Features
#### Custom Representations
You can customize the visual representation of the trajectory by passing rep_args during initialization. For example:

```Python
rep_args = [
    {"type": "cartoon", "params": {"color": "blue"}},
    {"type": "ball+stick", "params": {"color": "red"}}
]

image_generator = TrajectoryImageGenerator(
    traj=smoothed_traj, 
    image_folder='output_frames', 
    indices=range(0, len(smoothed_traj.trajectory), 10), 
    num_views=5, 
    rep_args=rep_args
)
```

### Atom Selection
Specify which atoms to visualize using the selection parameter. For example:

```Python
image_generator = TrajectoryImageGenerator(
    traj=smoothed_traj, 
    image_folder='output_frames', 
    indices=range(0, len(smoothed_traj.trajectory), 10), 
    num_views=5, 
    selection="protein or nucleic"
)
```

## Example Output
### Images
Generated images are saved in the specified image_folder with filenames like frame_000001.png, frame_000002.png, etc.

## Movie
The final movie is saved as an MP4 file, e.g., trajectory_movie.mp4.

## License
This project is licensed under the GNU Lesser General Public License v2.1. You may use, modify, and distribute this project under the terms of the LGPL v2.1. See the [LICENSE](LICENSE) file for details.

For more information about the LGPL v2.1, visit [https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html](https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html).

## Contributing
Contributions are welcome! If you encounter any issues or have feature requests, feel free to open an issue or submit a pull request.

## Acknowledgments
`MDAnalysis` for trajectory handling.

`NGLView` for visualization.

`MoviePy` for video creation.

## Author
This project was created by Emile de Bruyn.
