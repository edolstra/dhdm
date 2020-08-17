This is a tool to convert .dhdm files (HD Morphs) to tangent-space
displacement maps (encoded as 16-bit precision PNGs) that can be used
in Blender.

# Building / Running

If you're using the [Nix package manager](https://nixos.org/):

```console
# nix run github:edolstra/dhdm -- flags...
```

Otherwise, make sure you have the necessary dependencies (OpenSubdiv,
Boost, libpng, OpenGL, GLEW, GLFW and GLM) and build as follows:

```console
# cd src
# make -j
```

# Generating displacement maps

The general syntax is:

```console
# dhdm morphs-to-displacement <prefix> <base-mesh>.dsf <uv-map>.dsf [<morph>.[dsf|dhdm][=strength]]...
```

For example:

```console
# dhdm morphs-to-displacement \
    displacement \
    'data/daz 3d/genesis 8/female/genesis8female.dsf' \
    'data/daz 3d/genesis 8/female/uv sets/daz 3d/base/base female.dsf' \
    'fhmolympia8_hdlv4.dhdm' \
    'fbmolympia8_hdlv4.dhdm' \
    'vas forehead basic.dsf'=0.5
```

generates files named `displacement-{0,1,2,...}.png` that encode the
displacements given by the HD morphs 'fhmolympia8_hdlv4.dhdm',
'fbmolympia8_hdlv4.dhdm' and 'vas forehead basic.dsf', where the first
two have a strength of 100% and the last one has a strength of
50%. Each file corresponds to a UV tile; specifically,
`<prefix>-<N>.png` contains the UV faces with *u*-coordinates between
N and N+1. (All *v*-coordinates are assumed to be between 0 and 1.)
For the Genesis model this means: 0 = face; 1 = torso; 2 = legs; 3 =
arms; 4 = teeth; 5 = irises/pupils; 6 = cornea. No images are written
that don't contain any displaced geometry.

In addition to `.dhdm` HD morphs, you can also specify `.dsf`
morphs. However, it's generally better to use shape keys for such
"level 0" morphs. If you specify a .dsf morph and a .dhdm file with
the same base name also exists, they will both be used with the same
strength.

# Encoding

Displacement maps are stored as PNG files with 16 bits of precision
per channel. Each channel stores a displacement between -0.5 and 0.5
mapped to the range 0-1. (This means that with the default strengths,
the map can only store displacements up to 0.5 cm, but usually that's
sufficient.) The channels are as expected by Blender's Vector
Displacement node:

* Green: Displacement along the normal vector.
* Red: Displacement along the tangent vector.
* Blue: Displacement along the bitangent vector.

# Blender node setup

Connect the displacement map to a **Vector Displacement** node, and
connect its output to the Material Output's Displacement input, like
this:

![Node setup](/doc/node-setup.png)

Make sure of the following:

* Set the texture *Color Space* to **Non-Color**.

* Set the Vector Displacement node's *Type* to **Tangent Space**.

* Set the Vector Displacement node's *Midlevel* to **0.5**.

* Set the Vector Displacement node's *Scale* to whatever scale you
  used when you imported the model into Blender, e.g. when using 1
  Blender Unit (BU) = 1 cm, set it to 1.0; when using 1 BU = 1 m, set
  it to 0.01. Adjust as appropriate to increase or decrease the
  displacement effect.

* Make sure you're using the same UV map used to generate the
  displacement maps.

* The above setup assumes that no material spans multiple UV
  tiles. But if you have e.g. a single material, you'll need to mess
  around with Mapping nodes etc.

* Under the material's *Settings*, set *Displacement* to *Displacement
  and Bump* or *Displacement Only* to get true displacement.

* Add a subdivision modifier with either a sufficient number of
  subdivision levels to get the necessary geometry for displacement to
  work on (e.g. 3 or 4), or set it to "Adaptive".

* **Important**: Don't use normal maps in conjunction with vector
  displacement unless you know what you're doing! At least as of
  Blender 2.83, the "Normal Map" node does not appear to take
  displacements into account, so the resulting normals will be as if
  the displacements are not there (even though the geometry *is*
  displaced). This has the effect of largely obscuring the effects of
  (small) displacements, depending on the lighting conditions.

  However, you *can* use nodes like Bump that have a Normal input,
  since this allows them to use the displaced normal. So you can use
  bump maps on top of vector displacement maps as follows:

  ![Bump mapping node setup](/doc/bump-setup.png)
