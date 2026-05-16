# Metal Shaders for Immersed Boundary Flow Visualization

This folder contains Metal shaders for GPU-accelerated visualization of 2D CFD solutions produced by the MOVERS-N Euler / Immersed Boundary solvers.

## Files

- **FlowVisualization.metal** – Vertex and fragment shaders that map a 2D scalar field (density, pressure, Mach number, or velocity magnitude) to a colormap and optionally draw the immersed boundary in a solid color.
- **ContourOverlay.metal** – Fragment shader that draws iso-contour lines on the same scalar/mask layout; use with **flowVertex** and blend over the flow visualization.

## Data Layout

The solvers and `convert_to_vtk` output structured 2D data with:

- **Scalars:** `rho`, `pressure`, `mach`, velocity magnitude  
- **Vectors:** `u`, `v`  
- **Mask:** fluid vs solid (e.g. `Fluid_Region` in VTK, or `rho > 0`)

To use the shaders:

1. **Scalar texture** (texture 0): One channel (R) per grid point. Fill from your `.dat` or VTK (e.g. Mach, density, or pressure). Normalization is done in the shader using `valueMin` / `valueMax`.
2. **Mask texture** (texture 1): Same grid; R = 1 for fluid, 0 for solid/immersed. If you don’t use a mask, bind a 1×1 texture with value 1.0.
3. **FlowUniforms** (buffer 0): See below.

## FlowUniforms

| Field           | Type  | Description |
|----------------|-------|-------------|
| `valueMin`     | float | Min scalar value for colormap (e.g. min Mach in domain) |
| `valueMax`     | float | Max scalar value (e.g. max Mach) |
| `maskThreshold`| float | If no mask texture: treat as fluid when `scalar > maskThreshold` (e.g. 0 for density) |
| `colormapMode` | int   | 0 = jet, 1 = viridis, 2 = coolwarm, 3 = grayscale |
| `showBoundary` | int   | 1 = draw solid/immersed region with boundary color |
| `boundaryR/G/B`| float | Color for immersed boundary (0–1) |

## Entry Points

### FlowVisualization.metal

- **flowVertex** – Full-screen quad; buffers 0 = positions (`float2`), 1 = tex coords (`float2`).
- **flowFragment** – Main visualization: scalar texture + mask + FlowUniforms → RGBA.
- **flowFragmentVelocityMag** – Velocity-magnitude colormap only (same uniforms; texture 1 can be unused or repurposed).

### ContourOverlay.metal

- **contourFragment** – Iso-contour overlay. Use with **flowVertex** (same quad). Same textures: scalar (0), mask (1). Buffer 0 = **ContourUniforms**:

| Field          | Type  | Description |
|----------------|-------|-------------|
| `valueMin`     | float | Min scalar (same as FlowUniforms) |
| `valueMax`     | float | Max scalar |
| `maskThreshold`| float | Fluid when scalar > this |
| `numContours`  | int   | Number of contour levels (e.g. 16) |
| `lineWidth`    | float | Line thickness (normalized uv) |
| `lineR/G/B`    | float | Contour line color (0–1) |
| `invTexWidth`  | float | 1.0 / scalar texture width |
| `invTexHeight`| float | 1.0 / scalar texture height |

Render contour pass with alpha blending over the flow pass for combined visualization.

## Integration

- **macOS / iOS app:** Add `FlowVisualization.metal` to your Xcode target so it gets compiled into the default library. Create a `MTLRenderPipelineState` with the desired vertex/fragment functions and set textures and `FlowUniforms` each frame.
- **Standalone:** Compile with `metal -c FlowVisualization.metal -o FlowVisualization.air` and then `metallib FlowVisualization.air -o FlowVisualization.metallib` for loading at runtime.

Colormap and range are controlled via `FlowUniforms`; no shader recompile is needed to switch between density, pressure, or Mach number.
