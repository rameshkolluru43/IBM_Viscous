/*
 * FlowVisualization.metal
 * Implements specification in Shaders/README.md.
 * Metal shaders for 2D Immersed Boundary (IBM) CFD flow visualization.
 * Data layout: scalar texture (0), mask texture (1), FlowUniforms (buffer 0).
 */

#include <metal_stdlib>
using namespace metal;

// ---- Vertex output ----
struct VertexOut {
    float4 position [[position]];
    float2 texCoord;
};

// ---- FlowUniforms (buffer 0): matches README table ----
// valueMin, valueMax = colormap range; maskThreshold = fluid if scalar > this
// colormapMode: 0=jet, 1=viridis, 2=coolwarm, 3=grayscale
// showBoundary=1: draw solid region; boundaryR/G/B = boundary color (0-1)
struct FlowUniforms {
    float valueMin;
    float valueMax;
    float maskThreshold;
    int   colormapMode;
    int   showBoundary;
    float boundaryR;
    float boundaryG;
    float boundaryB;
};

// Full-screen quad vertex shader
vertex VertexOut flowVertex(uint vid [[vertex_id]],
                            constant float2* vertices [[buffer(0)]],
                            constant float2* texCoords [[buffer(1)]]) {
    VertexOut out;
    out.position = float4(vertices[vid], 0.0, 1.0);
    out.texCoord = texCoords[vid];
    return out;
}

// ---- Colormaps (s in [0,1]) ----
float4 colormapJet(float s) {
    s = clamp(s, 0.0, 1.0);
    float r = clamp(1.5 - 4.0 * fabs(s - 0.25), 0.0, 1.0);
    float g = clamp(1.5 - 4.0 * fabs(s - 0.5),  0.0, 1.0);
    float b = clamp(1.5 - 4.0 * fabs(s - 0.75), 0.0, 1.0);
    return float4(r, g, b, 1.0);
}

float4 colormapViridis(float s) {
    s = clamp(s, 0.0, 1.0);
    const float3 c0 = float3(0.267004, 0.004874, 0.329415);
    const float3 c1 = float3(0.282327, 0.140930, 0.457017);
    const float3 c2 = float3(0.127568, 0.566949, 0.550556);
    const float3 c3 = float3(0.369214, 0.788788, 0.383382);
    const float3 c4 = float3(0.993248, 0.906157, 0.143936);
    float t = s * 4.0;
    float3 c;
    if (t <= 1.0)      c = mix(c0, c1, t);
    else if (t <= 2.0) c = mix(c1, c2, t - 1.0);
    else if (t <= 3.0) c = mix(c2, c3, t - 2.0);
    else               c = mix(c3, c4, t - 3.0);
    return float4(c, 1.0);
}

float4 colormapCoolWarm(float s) {
    s = clamp(s, 0.0, 1.0);
    float r = 0.2298 + 0.7242 * s;
    float g = 0.2987 + 0.3917 * s - 0.8124 * s * s;
    float b = 0.7531 - 0.3984 * s;
    return float4(r, g, b, 1.0);
}

float4 colormapGrayscale(float s) {
    s = clamp(s, 0.0, 1.0);
    return float4(s, s, s, 1.0);
}

float4 applyColormap(float s, int mode) {
    switch (mode) {
        case 1: return colormapViridis(s);
        case 2: return colormapCoolWarm(s);
        case 3: return colormapGrayscale(s);
        default: return colormapJet(s);
    }
}

// Fragment shader: sample scalar texture, optional mask, apply colormap.
// scalarTex: single-channel (R) field (e.g. Mach, density, or pressure).
// maskTex:   fluid = 1, solid/immersed = 0. If not used, bind a 1x1 texture with value 1.
fragment float4 flowFragment(VertexOut in [[stage_in]],
                             texture2d<float, access::sample> scalarTex [[texture(0)]],
                             texture2d<float, access::sample> maskTex   [[texture(1)]],
                             constant FlowUniforms& uniforms [[buffer(0)]]) {
    constexpr sampler s(coord::normalized, filter::linear, address::clamp_to_edge);

    float2 uv = in.texCoord;
    float scalar = scalarTex.sample(s, uv).r;
    float mask = maskTex.sample(s, uv).r;
    bool isFluid = (mask > 0.5) || (scalar > uniforms.maskThreshold);

    if (!isFluid && uniforms.showBoundary != 0) {
        return float4(uniforms.boundaryR, uniforms.boundaryG, uniforms.boundaryB, 1.0);
    }

    float range = uniforms.valueMax - uniforms.valueMin;
    float t = (range > 1e-6) ? (scalar - uniforms.valueMin) / range : 0.0;
    return applyColormap(t, uniforms.colormapMode);
}

// ---- Alternative: velocity magnitude + vector overlay ----
// Second pipeline: scalar = velocity magnitude; optional second texture for u,v for glyphs
fragment float4 flowFragmentVelocityMag(VertexOut in [[stage_in]],
                                         texture2d<float, access::sample> scalarTex [[texture(0)]],
                                         texture2d<float, access::sample> uvTex     [[texture(1)]],
                                         constant FlowUniforms& uniforms [[buffer(0)]]) {
    constexpr sampler s(coord::normalized, filter::linear, address::clamp_to_edge);
    float2 uv = in.texCoord;
    float velMag = scalarTex.sample(s, uv).r;
    float s = (uniforms.valueMax > uniforms.valueMin)
        ? (velMag - uniforms.valueMin) / (uniforms.valueMax - uniforms.valueMin) : 0.0;
    return applyColormap(clamp(s, 0.0, 1.0), uniforms.colormapMode);
}
