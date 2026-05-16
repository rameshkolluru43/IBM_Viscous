/*
 * ContourOverlay.metal
 * Implements contour-line overlay for Shaders/README.md data layout.
 * Samples scalar texture (0), optional mask (1); draws iso-lines.
 * Use over the same full-screen quad as FlowVisualization (flowVertex).
 */

#include <metal_stdlib>
using namespace metal;

struct VertexOut {
    float4 position [[position]];
    float2 texCoord;
};

struct ContourUniforms {
    float valueMin;       // scalar range (same as FlowUniforms)
    float valueMax;
    float maskThreshold;
    int   numContours;    // number of contour levels (e.g. 16)
    float lineWidth;      // in normalized uv; typical 0.002–0.01
    float lineR, lineG, lineB;
    float invTexWidth;    // 1.0 / texture width (for gradient)
    float invTexHeight;   // 1.0 / texture height
};

// Fragment: draw contour lines at uniform scalar intervals
fragment float4 contourFragment(VertexOut in [[stage_in]],
                               texture2d<float, access::sample> scalarTex [[texture(0)]],
                               texture2d<float, access::sample> maskTex   [[texture(1)]],
                               constant ContourUniforms& uniforms [[buffer(0)]]) {
    constexpr sampler s(coord::normalized, filter::linear, address::clamp_to_edge);

    float2 uv = in.texCoord;
    float scalar = scalarTex.sample(s, uv).r;
    float mask = maskTex.sample(s, uv).r;
    bool isFluid = (mask > 0.5) || (scalar > uniforms.maskThreshold);

    if (!isFluid)
        return float4(0.0, 0.0, 0.0, 0.0); // transparent over solid

    float range = uniforms.valueMax - uniforms.valueMin;
    if (range < 1e-6)
        return float4(0.0, 0.0, 0.0, 0.0);

    float t = (scalar - uniforms.valueMin) / range;
    int n = max(1, uniforms.numContours);
    float level = floor(t * float(n) + 0.5) / float(n);

    float dx = uniforms.invTexWidth;
    float dy = uniforms.invTexHeight;
    float dtdx = (scalarTex.sample(s, uv + float2(dx, 0)).r - scalarTex.sample(s, uv - float2(dx, 0)).r) / (2.0 * range * max(dx, 1e-6));
    float dtdy = (scalarTex.sample(s, uv + float2(0, dy)).r - scalarTex.sample(s, uv - float2(0, dy)).r) / (2.0 * range * max(dy, 1e-6));
    float grad = length(float2(dtdx, dtdy)) + 0.001;
    float w = uniforms.lineWidth / grad;

    float dist = fabs(t - level);
    float alpha = 1.0 - smoothstep(0.0, w, dist);
    float3 lineColor = float3(uniforms.lineR, uniforms.lineG, uniforms.lineB);
    return float4(lineColor, alpha);
}
