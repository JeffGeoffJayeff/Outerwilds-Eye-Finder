"""
Plot points in spherical coordinates from npy files.
Takes polar and azimuth angles (elements 7 and 8) and plots them
with constant radius on a 3D sphere.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
from pathlib import Path
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import datashader as ds
import datashader.transfer_functions as tf

def load_npy_files(folder_path):
    """Load all npy files from a folder and return combined data."""
    npy_files = list(Path(folder_path).glob("*.npy"))
    
    if not npy_files:
        print(f"No .npy files found in {folder_path}")
        return None
    
    all_data = []
    for npy_file in npy_files:
        #print(f"Loading {npy_file.name}...")
        data = np.load(npy_file)
        #print(f"  Shape: {data.shape}, dtype: {data.dtype}")
        
        # Print field names if structured array
        #if data.dtype.names:
           # print(f"  Fields: {data.dtype.names}")
        
        all_data.append(data)
    
    # Combine all data
    try:
        combined_data = np.concatenate(all_data) if all_data else None
    except Exception as e:
        print(f"Error combining data: {e}")
        combined_data = None
    
    print(f"Total points loaded: {len(combined_data) if combined_data is not None else 0}")
    if combined_data is not None:
        print(f"Combined data shape: {combined_data.shape}")
    return combined_data

def analyze_probe_visits(data):
    """
    Analyze probe visits to different bodies and print statistics.
    
    Args:
        data: numpy structured array
    """
    if data is None or len(data) == 0:
        print("No data to analyze")
        return
    
    # Get all visit fields
    visit_fields = [
        'Sun Visits',
        'Sun Station Visits',
        'Ember Twin Visits',
        'Ash Twin Visits',
        'Timber Hearth Visits',
        'Attlerock Visits',
        'Brittle Hollow Visits',
        "Hollow's Lantern Visits",
        "Giant's Deep Visits",
        'Cannon Visits',
        'Dark Bramble Visits',
        'Interloper Visits',
        'White Hole Visits',
        'Stranger Visits',
        'Random Eye Visits',
        'Spacey Visits'
    ]
    
    # Print total visits for each body
    print("\n" + "=" * 70)
    print("TOTAL VISITS BY BODY (across all probes)")
    print("=" * 70)
    
    total_visits_by_body = {}
    for field in visit_fields:
        total = np.sum(data[field])
        total_visits_by_body[field] = total
        body_name = field.replace(' Visits', '')
        print(f"{body_name:30s}: {total:,} visits")
    
    total_all = sum(total_visits_by_body.values())
    print("-" * 70)
    print(f"{'TOTAL':30s}: {total_all:,} visits")
    
    # Find probes that visited multiple bodies
    print("\n" + "=" * 70)
    print("PROBES THAT VISITED MULTIPLE BODIES")
    print("=" * 70)
    
    multi_visit_probes = []
    
    for idx, probe in enumerate(data):
        # Count how many different bodies this probe visited
        bodies_visited = 0
        for field in visit_fields:
            if probe[field] > 0:
                bodies_visited += 1
        
        if bodies_visited > 2:
            multi_visit_probes.append((idx, probe, bodies_visited))
    
    print(f"\nFound {len(multi_visit_probes)} probes that visited multiple bodies\n")
    
    # Print launch conditions for multi-visit probes
    if multi_visit_probes:
        print("Launch conditions for probes visiting multiple bodies:")
        print("-" * 70)
        
        for probe_idx, probe, num_bodies in multi_visit_probes:
            launch_x = probe['Launch x']
            launch_y = probe['Launch y']
            launch_z = probe['Launch z']
            launch_vel = probe['Launch Velocity']
            
            # Calculate spherical coordinates of launch direction
            polar, azimuth, magnitude = cartesian_to_spherical(launch_x, launch_y, launch_z)
            
            # Count visits per body for this probe
            bodies_info = []
            for field in visit_fields:
                visits = probe[field]
                if visits > 0:
                    body_name = field.replace(' Visits', '')
                    bodies_info.append(f"{body_name} ({visits})")
            
            #print(f"\nProbe {probe_idx}:")
            #print(f"  Bodies visited: {', '.join(bodies_info)}")
            #print(f"  Launch direction (x,y,z): ({launch_x:.6f}, {launch_y:.6f}, {launch_z:.6f})")
            #print(f"  Launch direction (polar, azimuth): ({polar:.6f}, {azimuth:.6f})")
            #print(f"  Launch velocity: {launch_vel:.6f}")
    
    print("\n" + "=" * 70)

def spherical_to_cartesian(polar, azimuth, radius=1.0):
    """
    Convert spherical coordinates to Cartesian coordinates.
    
    Args:
        polar: Polar angle from north pole (0 to pi, or 0 to 180 degrees)
        azimuth: Azimuth angle in xy-plane (0 to 2*pi, or 0 to 360 degrees)
        radius: Radius from origin (default 1.0)
    
    Returns:
        x, y, z coordinates
    """
    x = radius * np.cos(polar) * np.sin(azimuth)
    y = radius * np.sin(polar) * np.sin(azimuth)
    z = radius * np.cos(azimuth)
    return x, y, z

def cartesian_to_spherical(x, y, z):
    """
    Convert Cartesian coordinates to spherical coordinates.
    
    Args:
        x, y, z: Cartesian components
    
    Returns:
        polar: Polar angle from north pole (0 to pi)
        azimuth: Azimuth angle in xy-plane (0 to 2*pi)
        radius: Distance from origin
    """
    radius = np.sqrt(x**2 + y**2 + z**2)
    
    # Handle zero radius
    radius = np.where(radius == 0, 1e-10, radius)
    
    polar = np.arccos(np.clip(z / radius, -1, 1))
    azimuth = np.arctan2(y, x)
    # Ensure azimuth is in [0, 2*pi]
    azimuth = np.where(azimuth < 0, azimuth + 2*np.pi, azimuth)
    
    return polar, azimuth, radius

def plot_spherical_heatmap(data, polar_field='Eye Shell Polar', azimuth_field='Eye Shell Azimuth', radius=1.0, bins=50):
    """
    Create a heatmap visualization on a 3D sphere using binned density.
    Much more efficient than scatter plot for large datasets.
    
    Args:
        data: numpy structured array
        polar_field: Field name for polar angle
        azimuth_field: Field name for azimuth angle
        radius: Radius of the sphere (default 1.0)
        bins: Number of bins for polar and azimuth (higher = more detail, slower)
    """
    if data is None or len(data) == 0:
        print("No data to plot")
        return
    
    print(f"Creating heatmap with {bins}x{bins} resolution")
    
    # Extract fields
    polar = np.array(data[polar_field], dtype=float)
    azimuth = np.array(data[azimuth_field], dtype=float)
    
    # Filter out NaN values
    valid_mask = np.isfinite(polar) & np.isfinite(azimuth)
    polar = polar[valid_mask]
    azimuth = azimuth[valid_mask]
    
    num_discarded = len(data) - len(polar)
    if num_discarded > 0:
        print(f"Discarded {num_discarded} points with NaN values")
    print(f"Using {len(polar)} valid points")
    
    if len(polar) == 0:
        print("No valid data points after filtering NaN values")
        return
    
    print(f"Polar range: [{polar.min():.4f}, {polar.max():.4f}]")
    print(f"Azimuth range: [{azimuth.min():.4f}, {azimuth.max():.4f}]")
    
    # Create 2D histogram (density map)
    heatmap, polar_edges, azimuth_edges = np.histogram2d(
        polar, azimuth, 
        bins=[bins, bins],
        range=[[polar.min(), polar.max()], [azimuth.min(), azimuth.max()]]
    )
    
    # Normalize heatmap for coloring
    heatmap_normalized = heatmap / heatmap.max()
    
    # Create grid for sphere surface
    polar_bins = np.linspace(polar.min(), polar.max(), bins)
    azimuth_bins = np.linspace(azimuth.min(), azimuth.max(), bins)
    
    polar_mesh, azimuth_mesh = np.meshgrid(polar_bins, azimuth_bins, indexing='ij')
    
    # Convert spherical to Cartesian for sphere surface
    x = radius * np.sin(polar_mesh) * np.cos(azimuth_mesh)
    y = radius * np.sin(polar_mesh) * np.sin(azimuth_mesh)
    z = radius * np.cos(polar_mesh)
    
    # Create 3D plot
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot surface colored by density
    surf = ax.plot_surface(x, y, z, facecolors=plt.cm.hot(heatmap_normalized),
                          shade=False, rstride=1, cstride=1)
    
    # Add colorbar
    sm = plt.cm.ScalarMappable(cmap=plt.cm.hot, 
                               norm=plt.Normalize(vmin=0, vmax=heatmap.max()))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.1, shrink=0.8)
    cbar.set_label('Point Density', rotation=270, labelpad=20)
    
    # Labels and formatting
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(f'Spherical Heatmap ({len(polar)} valid points, {bins}x{bins} bins)')
    
    # Set equal aspect ratio
    ax.set_box_aspect([1, 1, 1])
    
    plt.tight_layout()
    return fig, ax

def plot_launch_directions_heatmap(data, x_field='Launch x', y_field='Launch y', z_field='Launch z', radius=1.0, bins=50):
    """
    Create a heatmap of launch directions (Cartesian unit vectors) on a sphere.
    
    Args:
        data: numpy structured array
        x_field: Field name for x component
        y_field: Field name for y component
        z_field: Field name for z component
        radius: Radius of the sphere (default 1.0)
        bins: Number of bins for polar and azimuth (higher = more detail, slower)
    """
    if data is None or len(data) == 0:
        print("No data to plot")
        return
    
    print(f"Creating launch direction heatmap with {bins}x{bins} resolution")
    
    # Extract Cartesian components
    x = np.array(data[x_field], dtype=float)
    y = np.array(data[y_field], dtype=float)
    z = np.array(data[z_field], dtype=float)
    
    # Filter out NaN values
    valid_mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    x = x[valid_mask]
    y = y[valid_mask]
    z = z[valid_mask]
    
    num_discarded = len(data) - len(x)
    if num_discarded > 0:
        print(f"Discarded {num_discarded} points with NaN values")
    print(f"Using {len(x)} valid points")
    
    if len(x) == 0:
        print("No valid data points after filtering NaN values")
        return
    
    # Convert to spherical coordinates
    polar, azimuth, vector_magnitude = cartesian_to_spherical(x, y, z)
    
    print(f"Polar range: [{polar.min():.4f}, {polar.max():.4f}]")
    print(f"Azimuth range: [{azimuth.min():.4f}, {azimuth.max():.4f}]")
    print(f"Vector magnitude range: [{vector_magnitude.min():.4f}, {vector_magnitude.max():.4f}]")
    
    # Create 2D histogram (density map)
    heatmap, polar_edges, azimuth_edges = np.histogram2d(
        polar, azimuth, 
        bins=[bins, bins],
        range=[[0, np.pi], [0, 2*np.pi]]
    )
    
    # Normalize heatmap for coloring
    heatmap_normalized = heatmap / heatmap.max()
    
    # Create grid for sphere surface covering full sphere
    polar_bins = np.linspace(0, np.pi, bins)
    azimuth_bins = np.linspace(0, 2*np.pi, bins)
    
    polar_mesh, azimuth_mesh = np.meshgrid(polar_bins, azimuth_bins, indexing='ij')
    
    # Convert spherical to Cartesian for sphere surface
    x_sphere = radius * np.sin(polar_mesh) * np.cos(azimuth_mesh)
    y_sphere = radius * np.sin(polar_mesh) * np.sin(azimuth_mesh)
    z_sphere = radius * np.cos(polar_mesh)
    
    # Create 3D plot
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot surface colored by density
    surf = ax.plot_surface(x_sphere, y_sphere, z_sphere, facecolors=plt.cm.viridis(heatmap_normalized),
                          shade=False, rstride=1, cstride=1)
    
    # Add colorbar
    sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, 
                               norm=plt.Normalize(vmin=0, vmax=heatmap.max()))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.1, shrink=0.8)
    cbar.set_label('Point Density', rotation=270, labelpad=20)
    
    # Labels and formatting
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(f'Launch Directions Heatmap ({len(x)} valid points, {bins}x{bins} bins)')
    
    # Set equal aspect ratio
    ax.set_box_aspect([1, 1, 1])
    
    plt.tight_layout()
    return fig, ax

def plot_spherical_points(data, polar_field='Eye Shell Polar', azimuth_field='Eye Shell Azimuth', radius=1.0, range=100000):
    # Get Cartesian mesh grid
    sun_radius = 1000
    spherephi, spheretheta = np.mgrid[0.0:np.pi:20j, 0.0:2.0 * np.pi:20j] #Change the 20j to somethingelsej if you want different resolution on the sphere
    spherex = sun_radius*np.sin(spherephi) * np.cos(spheretheta) + 3421.723
    spherey = sun_radius*np.sin(spherephi) * np.sin(spheretheta) - 16097.95
    spherez = sun_radius*np.cos(spherephi)

    """
    Plot points in spherical coordinates on a 3D sphere.
    
    Args:
        data: numpy structured array
        polar_field: Field name for polar angle
        azimuth_field: Field name for azimuth angle
        radius: Radius of the sphere (default 1.0)
    """
    if data is None or len(data) == 0:
        print("No data to plot")
        return
    
    # Extract fields
    polar = np.array(data[polar_field], dtype=float)
    azimuth = np.array(data[azimuth_field], dtype=float)
    t = np.array(data['Eye Shell Time'], dtype=float)  # Assuming there's a time field for color mapping

    # Filter out NaN values in all relevant fields together
    valid_mask = np.isfinite(polar) & np.isfinite(azimuth) & np.isfinite(t)
    polar = polar[valid_mask]
    azimuth = azimuth[valid_mask]
    t = t[valid_mask]
    
    num_discarded = len(data) - len(polar)
    if num_discarded > 0:
        print(f"Discarded {num_discarded} points with NaN values")
    print(f"Using {len(polar)} valid points")
    
    if len(polar) == 0:
        print("No valid data points after filtering NaN values")
        return
    
    # Convert spherical to Cartesian coordinates
    x, y, z = spherical_to_cartesian(polar, azimuth, radius=radius)
    limit = range
    arr = np.zeros(len(x), dtype=bool)
    #arr[:limit] = True
    np.random.shuffle(arr)
    arr[:] = True
    df = pd.DataFrame({ #Take random samples 
        'x': x[arr],
        'y': y[arr],
        'z': z[arr],
        't': t[arr]
    })
    print(np.sum(np.isnan(df['x'])))
    print(np.sum(np.isnan(df['y'])))
    print(np.sum(np.isnan(df['z'])))
    print(np.sum(np.isnan(df['t'])))
    fig = px.scatter_3d(df, x='x', y='y', z='z', color='t')
    fig.update_scenes(aspectmode='cube')  # Making the axes be a cube
    fig.update_traces(marker_size=2)
    fig.add_surface(x=spherex, y=spherey, z=spherez, opacity=1.0,showscale=False)
    fig.show()

def plot_distance(data, polar_field='Eye Shell Polar', azimuth_field='Eye Shell Azimuth', radius='Final Sun Distance', range=100000):
    # Get Cartesian mesh grid for Giant's Deep
    GD_radius = 1000
    sun_radius = 2000
    spherephi, spheretheta = np.mgrid[0.0:np.pi:20j, 0.0:2.0 * np.pi:20j] #Change the 20j to somethingelsej if you want different resolution on the sphere
    GDx = GD_radius*np.sin(spherephi) * np.cos(spheretheta) + 3421.723
    GDy = GD_radius*np.sin(spherephi) * np.sin(spheretheta) - 16097.95
    GDz = GD_radius*np.cos(spherephi)

    sunx = sun_radius*np.sin(spherephi) * np.cos(spheretheta)
    suny = sun_radius*np.sin(spherephi) * np.sin(spheretheta)
    sunz = sun_radius*np.cos(spherephi)
    """
    Plot points in spherical coordinates on a 3D sphere.
    
    Args:
        data: numpy structured array
        polar_field: Field name for polar angle
        azimuth_field: Field name for azimuth angle
        radius: Radius of the sphere (default 1.0)
    """
    if data is None or len(data) == 0:
        print("No data to plot")
        return
    
    # Extract fields
    polar = np.array(data[polar_field], dtype=float)
    azimuth = np.array(data[azimuth_field], dtype=float)
    radius = np.array(data[radius], dtype=float)
    t = np.array(data['Eye Shell Time'], dtype=float)  # Assuming there's a time field for color mapping
    
    # Convert spherical to Cartesian coordinates
    x, y, z = spherical_to_cartesian(polar, azimuth, radius=radius)
    limit = range
    arr = np.zeros(len(x), dtype=bool)
    arr[:limit] = True
    np.random.shuffle(arr)
    #arr[:] = True
    df = pd.DataFrame({ #Take random samples 
        'x': x[arr],
        'y': y[arr],
        'z': z[arr],
        't': t[arr]
    })
    print(np.sum(np.isnan(df['x'])))
    print(np.sum(np.isnan(df['y'])))
    print(np.sum(np.isnan(df['z'])))
    print(np.sum(np.isnan(df['t'])))
    fig = px.scatter_3d(df, x='x', y='y', z='z', color='t')
    fig.update_scenes(aspectmode='cube')  # Making the axes be a cube
    fig.update_traces(marker_size=2)
    fig.add_surface(x=GDx, y=GDy, z=GDz, opacity=1.0,showscale=False)
    fig.add_surface(x=sunx, y=suny, z=sunz, opacity=1.0,showscale=False)
    fig.show()


def plot_launch_directions(data, x_field='Launch x', y_field='Launch y', z_field='Launch z', 
                           body_hit_colors=None,limit=100000):
    """
    Plot launch directions in 3D Cartesian coordinates.
    NaN points are colored discretely by 'Body Hit' index, valid points by 'Eye Shell Time' (continuous).
    
    Args:
        data: numpy structured array
        x_field: Field name for x component
        y_field: Field name for y component
        z_field: Field name for z component
        body_hit_colors: Dictionary mapping Body Hit index to color name. 
                        Example: {0: 'red', 1: 'blue', 2: 'green', 3: 'yellow'}
                        If None, uses default Plotly colors.
    """
    if data is None or len(data) == 0:
        print("No data to plot")
        return
    range = [-1,1]
    # Extract Cartesian components
    x = np.array(data[x_field], dtype=float)
    y = np.array(data[y_field], dtype=float)
    z = np.array(data[z_field], dtype=float)
    t = np.array(data['Eye Shell Time'], dtype=float)
    body = np.array(data['Body Hit'], dtype=float)

    arr = np.zeros(len(x), dtype=bool)
    arr[:limit] = True
    np.random.shuffle(arr)
    #arr[:] = True
    x = x[arr]
    y = y[arr]
    z = z[arr]
    t = t[arr]
    body = body[arr]

    print(f"Total points: {len(x)} \n Hit Something {np.sum(np.isfinite(body))} \n Hit Nothing {np.sum(np.isnan(body) & np.isnan(t))} \n Hit Eye Shell {np.sum(np.isfinite(t))}")

    # Create mask for valid points
    valid_mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    
    num_nan = np.sum(np.isnan(t))
    num_valid = np.sum(np.isfinite(t))
    if num_nan > 0:
        print(f"Found {num_nan} points with NaN values (will color by Body Hit)")
    print(f"Using {num_valid} valid points (colored by Eye Shell Time)")
    
    # Create figure
    fig = go.Figure()

    # Add valid points (colored by Eye Shell Time - continuous)
    if num_valid > 0:
        x_valid = x[np.isfinite(t)]
        y_valid = y[np.isfinite(t)]
        z_valid = z[np.isfinite(t)]
        t_valid = t[np.isfinite(t)]
        
        fig.add_trace(go.Scatter3d(
            x=x_valid, y=y_valid, z=z_valid,
            mode='markers',
            marker=dict(
                size=2,
                color=t_valid,
                colorscale='Viridis',
                showscale=True,
                colorbar=dict(title="Eye Shell Time", x=1.1),
            ),
            name='Valid Points',
            opacity=1
        ))
    # Endless Orbit:
    mask = np.isnan(t) & np.isnan(body)
    x_orbit = x[mask]
    y_orbit = y[mask]
    z_orbit = z[mask]
    fig.add_trace(go.Scatter3d(
        x=x_orbit,
        y=y_orbit,
        z=z_orbit,
        mode='markers',
        name="Orbiter",
        marker=dict(size=2, color='black'),
        legendgroup='orbit',
        showlegend=True,
    ))
    
    # Hit Bodies
    if num_nan > 0:
        x_nan = x[np.isfinite(body)]
        y_nan = y[np.isfinite(body)]
        z_nan = z[np.isfinite(body)]
        body_hit = body[np.isfinite(body)]
        
        # Map Body Hit indices to color names
        if body_hit_colors is None:
            # Default color sequence
            color_palette = ['red', 'blue', 'green', 'orange', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan']
            # Get unique non-NaN values
            unique_body_hits = np.unique(body_hit[~np.isnan(body_hit)])
            body_hit_colors = {int(i): color_palette[idx % len(color_palette)] 
                              for idx, i in enumerate(unique_body_hits)}
        
        # Group NaN points by Body Hit index for better legend
        df = pd.DataFrame({
            'x': x_nan, 
            'y': y_nan,
            'z': z_nan,
            'body': body_hit
        })
        for body, group in df.groupby('body'):
            fig.add_trace(go.Scatter3d(
                x=group['x'],
                y=group['y'],
                z=group['z'],
                mode='markers',
                name=f"hit {body}",
                marker=dict(size=2, color=body_hit_colors[body],opacity=0.5),
                legendgroup=body,   # ties it to the same legend entry
                showlegend=True,
            ))
    # TODO: Add plotting for when something neither hits a body nor reaches the eye shell, which is common due to the probe getting stuck in an orbit
    fig.update_layout(scene=dict(aspectmode='cube',xaxis=dict(range=range),yaxis=dict(range=range),zaxis=dict(range=range)))
    fig.show()

def plot_2D_crash(data,polar_field='Eye Shell Polar', azimuth_field='Eye Shell Azimuth'):
    canvas = ds.Canvas(plot_width=800*2, plot_height=400*2,x_range=(-np.pi,np.pi),y_range=(0,np.pi),x_axis_type='linear',y_axis_type='linear')
    polar = np.array(data[polar_field], dtype=float)
    azimuth = np.array(data[azimuth_field], dtype=float)
    t = np.array(data['Eye Shell Time'], dtype=float)  # Assuming there's a time field for color mapping

    # Filter out NaN values in all relevant fields together
    valid_mask = np.isfinite(polar) & np.isfinite(azimuth) & np.isfinite(t)
    polar = polar[valid_mask]
    azimuth = azimuth[valid_mask]
    t = t[valid_mask]
    
    num_discarded = len(data) - len(polar)
    if num_discarded > 0:
        print(f"Discarded {num_discarded} points with NaN values")
    print(f"Using {len(polar)} valid points")
    
    if len(polar) == 0:
        print("No valid data points after filtering NaN values")
        return
    
    # Convert spherical to Cartesian coordinates
    limit = 500000
    df = pd.DataFrame({
        'polar':polar,
        'azi': azimuth,
        't': t
    })

    agg = canvas.points(df, 'polar', 'azi', agg=ds.count())
    ds.tf.set_background(ds.tf.shade(agg, cmap=plt.cm.viridis), "white").to_pil().save("2D_crash.png")


def main():

    """Main execution function."""
    outputs_folder = Path("OutputsWithDistance")
    
    if not outputs_folder.exists():
        print(f"Outputs folder not found at {outputs_folder.absolute()}")
        print("Creating a demo using existing npy files instead...")
        
        # Try with existing npy files in current directory
        current_folder = Path(".")
        data = load_npy_files(current_folder)
    else:
        data = load_npy_files(outputs_folder)
    
    if data is not None:
        # Analyze visits
        analyze_probe_visits(data)
        
        # Auto-select bin resolution based on dataset size
        num_points = len(data)
        bins = min(100, max(30, int(np.sqrt(num_points / 100))))
        print(f"\nAuto-selecting {bins} bins for {num_points} points\n")
        if False:
        # Plot 1: Eye Shell position heatmap
            print("=" * 60)
            print("PLOT 1: Eye Shell Position Heatmap")
            print("=" * 60)
            fig1, ax1 = plot_spherical_heatmap(data, 
                                            polar_field='Eye Shell Polar', 
                                            azimuth_field='Eye Shell Azimuth', 
                                            radius=1.0,
                                            bins=bins)
            
            # Save figure
            output_file1 = "spherical_heatmap_eye_shell.png"
            plt.savefig(output_file1, dpi=150, bbox_inches='tight')
            print(f"Plot saved to {output_file1}\n")
        if False:
            # Plot 2: Launch directions heatmap
            print("=" * 60)
            print("PLOT 2: Launch Directions Heatmap")
            print("=" * 60)
            fig2, ax2 = plot_launch_directions_heatmap(data,
                                                    x_field='Launch x',
                                                    y_field='Launch y',
                                                    z_field='Launch z',
                                                    radius=1.0,
                                                    bins=bins)
            
            # Save figure
            output_file2 = "spherical_heatmap_launch_directions.png"
            plt.savefig(output_file2, dpi=150, bbox_inches='tight')
            print(f"Plot saved to {output_file2}\n")
            
            plt.show()

        plot_spherical_points(data, polar_field='Eye Shell Polar', azimuth_field='Eye Shell Azimuth', radius=286500,range=400000)
        
        # Define custom colors for Body Hit indices (NaN points only)
        colormap = {
            0:"#FFDF22",
            1:"#9D00FF",
            2:"#f05f44",
            3:"#e19d49",
            4: "#5a8240",
            5: "#8e7263",
            6: "#647297",
            7: "#d86f23",
            8: "#31a174",
            9:"#676e4c",
            10: "#55503a",
            11: "#349fb9",
            12: "#D3D3D3",
            13: "#22a185",
            14: "#8673A1",
            15: "#00FF6A"
        }
        plot_launch_directions(data, x_field='Launch x', y_field='Launch y', z_field='Launch z',body_hit_colors=colormap,limit=400000)
        plot_distance(data, polar_field='Final Polar', azimuth_field='Final Azimuth', radius='Final Radius', range=400000)
        plot_2D_crash(data, polar_field='Eye Shell Polar', azimuth_field='Eye Shell Azimuth')
    else:
        print("No data was loaded.")

if __name__ == "__main__":
    main()
