# ##### BEGIN GPL LICENSE BLOCK #####
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software Foundation,
#  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# ##### END GPL LICENSE BLOCK #####

import bpy
import bmesh
import sys
import time
import os
from bpy.types import PropertyGroup, Operator, Panel, AddonPreferences
from mathutils import Vector, Matrix
from bpy_extras.io_utils import ImportHelper
from collections import Counter
from bpy.props import (
        BoolProperty,
        PointerProperty,
        FloatProperty,
        StringProperty,
        EnumProperty,
        IntProperty,
        FloatVectorProperty,
        )

bl_info = {
    "name": "kei2m",
    "author": "Kjell Emanuelsson",
    "category": "Import-Export",
    "version": (1, 3, 0, 6),
    "blender": (2, 80, 0),
    "location": "Viewport / N-Panel / kei2m",
    "warning": "",
    "description": "Image(s) To Mesh Generator",
    "doc_url": "https://ke-code.xyz",
}


kei2m_version = 1.306


def is_bversion(req_ver):
    bv = int("".join([str(i) for i in bpy.app.version]))
    if bv < req_ver:
        return False
    return True


def make_entry(h, width, w, pixels):
    idx = (h * width) + w
    px_index = idx * 4
    return pixels[px_index:px_index + 4]


def alpha_check(images, rgb=False, c2m=False):
    has_alpha = True
    for img in images:
        if img.depth != 32 and not rgb and not c2m:
            has_alpha = False
    return has_alpha


def reduce_colors(pixelmap, threshold=0.51, cap=None):
    colors = [p[2] for p in pixelmap if p[2][3] == 1]
    color_groups = []
    common_colors = [c for c in Counter(colors).most_common()]
    if cap is not None:
        if len(common_colors) > cap:
            common_colors = common_colors[:cap]
    for c in common_colors:
        rmin = c[0][0] - threshold
        rmax = c[0][0] + threshold
        gmin = c[0][1] - threshold
        gmax = c[0][1] + threshold
        bmin = c[0][2] - threshold
        bmax = c[0][2] + threshold
        similar = [[c][0]]
        for oc in reversed(common_colors):
            if oc != c:
                if rmin < oc[0][0] < rmax and gmin < oc[0][1] < gmax and bmin < oc[0][2] < bmax:
                    similar.append(oc[0])
                    common_colors.remove(oc)
        color_groups.append(similar)
    for p in pixelmap:
        for i, g in enumerate(color_groups):
            if p[2] in g:
                p[2] = g[0][0]
    return pixelmap, [g[0][0] for g in color_groups]


class KeI2M(Operator):
    bl_idname = "ke.i2m"
    bl_label = "kei2m"
    bl_description = "Image(s) To Mesh Generator"
    bl_options = {'REGISTER', 'UNDO'}

    opacity: IntProperty(min=1, max=100, default=50, name="Opacity Tolerance", subtype="PERCENTAGE",
                         soft_min=50, soft_max=50,
                         description="100% = Will only allow completely opaque image pixels\n"
                                     "Do not use as slider - Use keyboard input!")
    workres: EnumProperty(
        items=[("64", "64x64", "", "", 1),
               ("128", "128x128", "", "", 2),
               ("256", "256x256", "", "", 3),
               ("512", "512x512", "", "", 4),
               ("1024", "1024x1024", "", "", 5),
               ("IMAGE", "Image Size", "", "", 6),
               ],
        name="Work Resolution", default="128",
        description="The internal res. used for conversion.\nUse as low value as possible,\n"
                    "1k or more can be VERY slow. (The complexity of the alpha will also matter)\n")

    custom_workres: IntProperty(min=0, max=65536, default=0, name="Custom Resolution",
                                soft_min=0, soft_max=0,
                                description="Non-zero value will override Work Resolution X & Y sizes")

    screw_flip: BoolProperty(default=False, name="Flip Screw Source-Side",
                             description="Flips the side used in a Screw geo conversion")

    screw_xcomp: IntProperty(min=-100, max=100, default=15, name="Stretch Comp.", subtype="PERCENTAGE",
                             description="Offset X projection scale to compensate for screw projection distortion\n"
                                         "Do not use as slider! Use keyboard input!")
    reduce: EnumProperty(
        items=[("REDUCED", "Reduced", "", "", 1),
               ("SIMPLE", "Simple", "", "", 2),
               ("DISSOLVE", "Dissolve", "", "", 3),
               ("NONE", "None", "", "", 4),
               ],
        name="Mesh Reduction", default="SIMPLE",
        description="Which type of mesh reduction to use:\n"
                    "Reduced: No Smoothing & high polycount\n"
                    "Simple: Low smoothing & high polycount\n"
                    "Dissolve: Low polycount & high smoothness & reduction, slow\n"
                    "None: Fast, but VERY high polycount. 1 workpixel = 1 face!\n")

    c2m_reduce: EnumProperty(
        items=[("DISSOLVE", "Dissolve", "", "", 1),
               ("NONE", "None", "", "", 2),
               ],
        name="Mesh Reduction", default="DISSOLVE",
        description="Which type of mesh reduction to use:\n"
                    "Dissolve: Low polycount, Slow\n"
                    "None: High polycount, Fast. 1 workpixel = 1 face!\n")

    shade_smooth: BoolProperty(default=False, name="Smooth Shading", description="Smooth Shading & Auto-Smooth")

    front_only: BoolProperty(default=False, name="Front Texture Only",
                             description="Only use Front Image for Texturing,\n"
                                         "even if using Right/Top Boolean.\n"
                                         "(So the Right/Top images can be just silhouettes)")

    c2threshold: FloatProperty(min=0, max=1, default=0.52, name="Color Threshold",
                               description="Tolerance for color separation / reduction (anti-aliasing removal)\n"
                                           "Do not use as slider! Use keyboard input!")

    c2m_smooth: IntProperty(min=0, max=100, default=100, name="Mesh Smoothing", subtype="PERCENTAGE",
                            soft_min=100, soft_max=100,
                            description="Smooths vertices below the pixel edge length threshold\n"
                                        "(Avoiding long straight edges)\n"
                                        "Do not use as slider! Use keyboard input!")

    qnd_mat: BoolProperty(default=False, name="QnD Materials",
                          description="Will add Quick-n-Dirty Roughness and Bump nodes in addition to color.\n"
                                      "Based on the color input.")

    apply: BoolProperty(default=False, name="Apply All Modifiers",
                        description="Applies the Screw (if used) & UV Projection Modifiers and removes the Empties)\n"
                                    "Set to off to leave modifier live (for live uv projection modeling etc)")

    apply_none: BoolProperty(default=False, name="Dont Apply Modifiers",
                             description="Boolean Geo Mode: leaves Solidify & Boolean Modifiers active\n"
                                         "For debugging & troubleshooting mostly")

    angle: FloatProperty(min=0, max=0.9, default=0.5, name="UV Angle Tolerance",
                         soft_min=0.5, soft_max=0.5,
                         description="0 to 0.9 Angle tolerance for axis-switching (from one axis to the other)")

    reset: BoolProperty(default=False, name="Reset", description="Reset i2m to default values",
                        options={"SKIP_SAVE", "HIDDEN"})

    batch: BoolProperty(default=False, name="Batch", description="Batch process all images in a folder with i2m",
                        options={"SKIP_SAVE", "HIDDEN"})

    vcolor: BoolProperty(default=False, name="Vertex Color",
                         description="aka 'Retro Pixel Graphics'\n"
                                     "Each Pixel makes up one face, "
                                     "using Vertex Color instead of the Image Texture.\n"
                                     "Calculated from the Image Texture")

    vcthreshold: FloatProperty(min=0, max=1, default=0, name="Color Threshold",
                               soft_min=0, soft_max=0,
                               description="Tolerance for color separation / reduction (anti-aliasing removal)\n"
                                           "0 = No limit (full rgb) in Vertex Color Mode (Also Faster)\n"
                                           "Sensitive: Increase by steps of 0.05 (Also, very slow!)\n"
                                           "Do not use as slider! Use keyboard input!")

    dilation: IntProperty(min=0, max=99, default=0, name="Expand Border",
                          soft_min=0, soft_max=0,
                          description="Using a simple brute force dilation on the alpha (in pixels, roughly),\n"
                                      "expanding it beyond the original image alpha borders. Zero to disable.\n"
                                      "Tip: Tweak Tolerance value 1st - only use EB if necessary")

    pixel_width: FloatProperty(min=0, max=1, default=0, name="Pixel Width", precision=5,
                               soft_min=0, soft_max=0,
                               description="Set pixel size in BU (meter)\n"
                                           "Mesh Width is calculated with Pixel Width * Work Res")

    width: FloatProperty(min=0, default=1, name="Mesh Width", precision=3,
                         soft_min=1, soft_max=1,
                         description="Set Custom width/size in BU (meter)")

    size: EnumProperty(
        items=[("AUTOFIT", "Autofit", "", "", 1),
               ("WIDTH", "Custom Width", "", "", 2),
               ("PIXEL", "Pixel Width", "", "", 3),
               ],
        name="Size", default="AUTOFIT",
        description="Autofit: Automatically fits all resolutions to 1 BU (meter)\n"
                    "Width: Set custom mesh size/width\n"
                    "Pixel: Calculate size based on set custom pixel size/width \n")

    coll = None
    wm = None
    t = 0
    tot = 0
    img_count = 0
    flip_left = False
    flip_back = False
    noz = False
    use_rgb = False
    rgb = (1, 1, 1)
    cmats = []
    has_alpha = True
    color_cap = 16
    geo = "PLANE"
    c2m = False

    @classmethod
    def poll(cls, context):
        return context.mode == "OBJECT"

    def draw(self, context):
        k = context.scene.kei2m
        c2m_mode = True if k.geo == "C2M" else False

        layout = self.layout
        layout.use_property_split = True

        if c2m_mode:
            row = layout.row().split(factor=0.33)
            row.enabled = False
            row.separator()
            row.label(text="Color 2 Material Mode")
        elif self.vcolor:
            row = layout.row().split(factor=0.33)
            row.separator()
            row.enabled = False
            row.label(text="Vertex Color")

        if self.has_alpha:
            layout.prop(self, "opacity")
            layout.prop(self, "dilation")
            layout.separator(factor=0.5)

        if self.custom_workres == 0:
            layout.prop(self, "workres")

        layout.prop(self, "custom_workres")
        # if self.geo != "BOOLEAN":
        layout.separator(factor=0.5)
        layout.prop(self, "size", expand=True)
        if self.size == "WIDTH":
            layout.prop(self, "width")
        elif self.size == "PIXEL":
            layout.prop(self, "pixel_width")
        layout.separator(factor=0.5)

        if not self.vcolor:
            layout.prop(self, "reduce", expand=True)
            layout.separator(factor=0.5)

        # mode specifics
        if c2m_mode:
            layout.prop(self, "c2threshold", expand=True)
            layout.prop(self, "c2m_reduce", expand=True)
            if self.c2m_reduce != "NONE":
                layout.prop(self, "c2m_smooth", expand=True)
            layout.separator(factor=0.5)
        elif k.geo == "SCREW":
            layout.prop(self, "screw_xcomp", toggle=True)
            layout.prop(self, "screw_flip", toggle=True)
            layout.separator(factor=0.5)
        elif k.geo == "BOOLEAN":
            layout.prop(self, "angle")
            layout.prop(self, "front_only", toggle=True)
            layout.separator(factor=0.5)

        layout.prop(self, "vcolor", toggle=True)
        if self.vcolor:
            layout.prop(self, "vcthreshold", expand=True)
        
        if not is_bversion(410):
            layout.prop(self, "shade_smooth", toggle=True)
        
        if not c2m_mode:
            if not self.vcolor:
                layout.prop(self, "qnd_mat", toggle=True)
            if not self.apply and k.geo == "BOOLEAN":
                layout.prop(self, "apply_none", toggle=True)
            if not self.apply_none:
                layout.prop(self, "apply", toggle=True)
        row = layout.row(align=False)
        row.alignment = "LEFT"
        row.operator("wm.operator_defaults", icon="FILE_REFRESH", text="Reset")
        layout.separator()

    def make_pixel_map(self, width, height, pixels, scl, use_rgb=False):
        pixel_map = []
        width_range = width
        start = 0
        if self.geo == "SCREW":
            if not self.screw_flip:
                width_range = int(width / 2)
            else:
                start = int(width / 2)

        tolerance = float(self.opacity / 100)
        # hard to find opc value that "feels" good here...
        opc = self.opacity * 0.5
        tol = float(opc / 100)
        rmin = self.rgb[0] - tol
        rmax = self.rgb[0] + tol
        gmin = self.rgb[1] - tol
        gmax = self.rgb[1] + tol
        bmin = self.rgb[2] - tol
        bmax = self.rgb[2] + tol

        if self.dilation != 0:
            # Dilate alpha border
            tot = len(pixels)
            pm_dilated = list(pixels)

            for w in range(start, width_range):
                for h in range(0, height):
                    rgba = make_entry(h, width, w, pixels)
                    if use_rgb:
                        if rmin < rgba[0] < rmax and gmin < rgba[1] < gmax and bmin < rgba[2] < bmax:
                            pass
                        else:
                            for dx in range(-self.dilation, self.dilation):
                                for dy in range(-self.dilation, self.dilation):
                                    i = (((h + dy) * width) + (w + dx)) * 4
                                    if 0 < i < tot:
                                        pm_dilated[i:i + 3] = (0, 0, 0)
                    else:
                        if rgba[3] >= tolerance:
                            for dx in range(-self.dilation, self.dilation):
                                for dy in range(-self.dilation, self.dilation):
                                    i = (((h + dy) * width) + (w + dx)) * 4
                                    if 0 < i < tot:
                                        pm_dilated[i + 3] = 1
            pixels = tuple(pm_dilated)

        # Apply alpha tolerance (trim outline)
        if use_rgb:
            for w in range(start, width_range):
                for h in range(0, height):
                    rgba = make_entry(h, width, w, pixels)
                    if rmin < rgba[0] < rmax and gmin < rgba[1] < gmax and bmin < rgba[2] < bmax:
                        pass
                    else:
                        pixel_map.append([w * scl, h * scl, rgba])
        else:
            for w in range(start, width_range):
                for h in range(0, height):
                    rgba = make_entry(h, width, w, pixels)
                    if rgba[3] >= tolerance:
                        pixel_map.append([w * scl, h * scl, rgba])
        # Limit colors
        if self.c2m:
            pixel_map, self.cmats = reduce_colors(pixelmap=pixel_map, threshold=self.c2threshold, cap=self.color_cap)
        elif self.vcolor and self.vcthreshold > 0:
            pixel_map, self.cmats = reduce_colors(pixelmap=pixel_map, threshold=self.vcthreshold, cap=None)

        return pixel_map

    def make_mesh_data(self, pixel_map, work_res, scl, name, axis="Front"):
        s = scl * 0.5
        w = (work_res * scl) * 0.5
        verts = []
        faces = []
        if axis == "Right":
            for i, v in enumerate(pixel_map):
                x, y = v[:2]
                x -= w
                verts.extend([[0, -s + x, s + y], [0, -s + x, -s + y], [0, s + x, -s + y], [0, s + x, s + y]])
                offset = i * 4
                faces.append([0 + offset, 1 + offset, 2 + offset, 3 + offset])
        elif axis == "Top":
            for i, v in enumerate(pixel_map):
                x, y = v[:2]
                x -= w
                y -= w
                verts.extend(
                    [[-s + x, s + y, w], [-s + x, -s + y, w], [s + x, -s + y, w], [s + x, s + y, w]])
                offset = i * 4
                faces.append([0 + offset, 1 + offset, 2 + offset, 3 + offset])
        else:  # Front
            for i, v in enumerate(pixel_map):
                x, y = v[:2]
                x -= w
                verts.extend([[-s + x, 0, s + y], [-s + x, 0, -s + y], [s + x, 0, -s + y], [s + x, 0, s + y]])
                offset = i * 4
                faces.append([0 + offset, 1 + offset, 2 + offset, 3 + offset])

        mesh = bpy.data.meshes.new(name)
        mesh.from_pydata(verts, [], faces)

        if self.screw_flip:
            mesh.flip_normals()

        if self.vcolor or self.c2m:
            vc = mesh.vertex_colors.new()

            for f, c in zip(mesh.polygons, pixel_map):
                if self.vcolor:
                    for loop in f.loop_indices:
                        vc.data[loop].color = c[2]
                if self.c2m:
                    for i, color in enumerate(self.cmats):
                        if c[2] == color:
                            f.material_index = i

        mesh.update()
        return mesh

    def cleanup(self, mesh, scl, axis="Front"):
        bm = bmesh.new()
        bm.from_mesh(mesh)
        bmesh.ops.remove_doubles(bm, verts=bm.verts, dist=scl * 0.25)
        inner_verts = [v for v in bm.verts if not v.is_boundary]

        if self.reduce == "DISSOLVE":
            scl_max = scl * 1.1
            smoothverts = []

            if self.c2m:
                bmesh.ops.dissolve_limit(bm, angle_limit=0.08728, verts=bm.verts, edges=bm.edges, delimit={"MATERIAL"})
                if self.c2m_smooth != 0:
                    smedges = []
                    for e in bm.edges:
                        if e.calc_length() < scl_max:
                            if len(e.verts[0].link_edges) < 3:
                                smoothverts.append(e.verts[0])
                            if len(e.verts[1].link_edges) < 3:
                                smoothverts.append(e.verts[1])
                            smedges.append(e)
                    smoothverts = list(set(smoothverts))
                    bmesh.ops.smooth_vert(bm, verts=smoothverts, factor=(self.c2m_smooth / 100) * 0.5,
                                          use_axis_x=True, use_axis_y=True, use_axis_z=True)
            else:
                for f in bm.faces:
                    es = []
                    for e in f.edges:
                        if e.is_boundary:
                            if e.calc_length() < scl_max:
                                es.extend(e.verts[:])
                    if len(es) == 4:
                        visited = set()
                        dup = {v for v in es if v in visited or (visited.add(v) or False)}
                        if dup:
                            smoothverts.extend(dup)
                smoothverts = list(set(smoothverts))
                if smoothverts:
                    bmesh.ops.dissolve_verts(bm, verts=smoothverts, use_face_split=True, use_boundary_tear=True)
                    bmesh.ops.unsubdivide(bm, verts=inner_verts, iterations=64)
                    bmesh.ops.dissolve_limit(bm, angle_limit=0.08727, verts=bm.verts, edges=bm.edges)

        elif self.reduce == "SIMPLE":
            bmesh.ops.unsubdivide(bm, verts=inner_verts, iterations=64)
            if not self.geo == "BOOLEAN":
                smoothverts = [v for v in bm.verts if v.is_boundary]
                bmesh.ops.smooth_vert(bm, verts=smoothverts, factor=1.0,
                                      use_axis_x=True, use_axis_y=True, use_axis_z=True)
            else:
                bmesh.ops.dissolve_limit(bm, angle_limit=0.08727, verts=bm.verts, edges=bm.edges)

        elif self.reduce == "REDUCED":
            bmesh.ops.unsubdivide(bm, verts=inner_verts, iterations=64)

        # ELSE: NO reduction - FULL PIXELATION

        # COMPENSATE SCALE OFFSET "FIX"
        c = scl * 0.5
        mtx = Matrix()
        if axis == "Front":
            bmesh.ops.translate(bm, vec=Vector((c, 0, c)), space=mtx, verts=bm.verts)
        elif axis == "Right":
            bmesh.ops.translate(bm, vec=Vector((0, c, c)), space=mtx, verts=bm.verts)
        elif axis == "Top":
            bmesh.ops.translate(bm, vec=Vector((c, c, 0)), space=mtx, verts=bm.verts)

        bm.to_mesh(mesh)
        bm.free()

    def make_scene_object(self, mesh, name):
        obj = bpy.data.objects.new(name, mesh)
        self.coll.objects.link(obj)
        if self.shade_smooth:
            obj.data.use_auto_smooth = True
        return obj

    def make_projector(self, p_name, w, axis="Front"):
        existing = bpy.data.objects.get(p_name)
        if existing:
            bpy.data.objects.remove(existing)
        projector = bpy.data.objects.new(p_name, None)
        self.coll.objects.link(projector)
        projector.hide_viewport = True
        projector.empty_display_type = 'SINGLE_ARROW'
        if axis == "Right":
            projector.rotation_euler = [-1.5707963, 3.1415926, -1.5707963]
        elif axis == "Top":
            projector.rotation_euler[2] = 1.5707963
        elif axis == "Back":
            projector.rotation_euler = [-1.5707963, 3.1415926, 0]
        elif axis == "Left":
            projector.rotation_euler = [1.5707963, 0, -1.5707963]
        elif axis == "Bottom":
            projector.rotation_euler = [-3.1415926, 0, 1.5707963]
        else:  # Front
            projector.rotation_euler[0] = 1.5707963
        projector.location[2] = w
        # Hack compensation
        if self.geo == "SCREW":
            val = (self.screw_xcomp / 1000) + w
            projector.scale = (val, w, w)
        else:
            projector.scale = (w, w, w)
        return projector

    def make_material(self, name, img):
        m = bpy.data.materials.get(name)
        if m:
            bpy.data.materials.remove(m)
        m = bpy.data.materials.new(name=name)
        m.use_nodes = True
        shader = m.node_tree.nodes["Material Output"].inputs[0].links[0].from_node
        # Color
        n_color = m.node_tree.nodes.new("ShaderNodeTexImage")
        m.node_tree.links.new(shader.inputs["Base Color"], n_color.outputs[0])
        n_color.image = img
        n_color.location = (-1100, 215)
        # Quick-and-Dirty Colormap-Texturing
        if self.qnd_mat:
            n_rbgmix = m.node_tree.nodes.new("ShaderNodeMixRGB")
            n_rbgmix.location = (-800, 100)
            m.node_tree.links.new(n_color.outputs[0], n_rbgmix.inputs[1])

            n_rbg2bw = m.node_tree.nodes.new("ShaderNodeRGBToBW")
            n_rbg2bw.location = (-600, 100)
            m.node_tree.links.new(n_rbgmix.outputs[0], n_rbg2bw.inputs[0])

            n_invert = m.node_tree.nodes.new("ShaderNodeInvert")
            n_invert.location = (-400, 25)
            m.node_tree.links.new(n_rbg2bw.outputs[0], n_invert.inputs[1])

            n_bump = m.node_tree.nodes.new("ShaderNodeBump")
            n_bump.inputs[0].default_value = 0.05
            n_bump.inputs[1].default_value = 0.05
            n_bump.location = (-400, -265)
            m.node_tree.links.new(n_rbg2bw.outputs[0], n_bump.inputs[2])
            m.node_tree.links.new(shader.inputs["Normal"], n_bump.outputs[0])

            n_mul = m.node_tree.nodes.new("ShaderNodeMath")
            n_mul.operation = "MULTIPLY"
            n_mul.location = (-200, 22)
            m.node_tree.links.new(n_invert.outputs[0], n_mul.inputs[1])
            m.node_tree.links.new(shader.inputs["Roughness"], n_mul.outputs[0])
        return m

    def sort_material_slots(self, m_axis):
        # Project "through" if axis is missing (front+back etc)
        idx = []
        vecs = []
        for i, a in enumerate(m_axis):
            if a == "Front":
                idx.append(i)
                vecs.append(Vector((0, 1, 0)))
                if "Back" not in m_axis:
                    self.flip_back = True
                    idx.append(i)
                    vecs.append(Vector((0, 1, 0)))
            elif a == "Right":
                idx.append(i)
                vecs.append(Vector((-1, 0, 0)))
                if "Left" not in m_axis:
                    self.flip_left = True
                    idx.append(i)
                    vecs.append(Vector((1, 0, 0)))
            elif a == "Top":
                idx.append(i)
                vecs.append(Vector((0, 0, -1)))
                if "Bottom" not in m_axis:
                    idx.append(i)
                    vecs.append(Vector((0, 0, 1)))
            elif a == "Back":
                idx.append(i)
                vecs.append(Vector((0, -1, 0)))
                if "Front" not in m_axis:
                    idx.append(i)
                    vecs.append(Vector((0, 1, 0)))
            elif a == "Left":
                idx.append(i)
                vecs.append(Vector((1, 0, 0)))
                if "Right" not in m_axis:
                    idx.append(i)
                    vecs.append(Vector((-1, 0, 0)))
            elif a == "Bottom":
                idx.append(i)
                vecs.append(Vector((0, 0, 1)))
                if "Top" not in m_axis:
                    idx.append(i)
                    vecs.append(Vector((0, 0, -1)))
        return idx, vecs

    def progress_update(self, context, txt, done):
        if done:
            tstat = round((time.time() - self.t), 5)
            t = '{:f}'.format(tstat).rstrip('0')
            if len(t) > 6:
                t = t[:6]
            msg = "\r{0}: [   COMPLETE   ] {1}s\r\n".format(txt, t)
            self.tot += float(t)
            self.wm.progress_update(99)
        else:
            msg = "\r{0}: [ Processing...]".format(txt)
            self.wm.progress_update(98)
        sys.stdout.write(msg)
        sys.stdout.flush()
        self.t = time.time()

    def execute(self, context):
        self.tot = 0

        if self.reset:
            self.opacity = 95
            self.workres = "128"
            self.geo = "PLANE"
            self.screw_flip = False
            self.screw_xcomp = 15
            self.reduce = "SIMPLE"
            self.shade_smooth = False
            self.front_only = False
            self.qnd_mat = False
            self.apply = False
            self.apply_none = False
            self.angle = 0.5
            self.reset = False
            self.custom_workres = 0
            self.vcolor = False
            return {"FINISHED"}

        k = context.scene.kei2m
        kap = context.preferences.addons["ke_i2m"].preferences

        self.color_cap = kap.cap

        self.use_rgb = kap.use_rgb
        if self.use_rgb:
            self.rgb = kap.user_rgb

        self.geo = k.geo

        # Forced Overrides
        if self.geo == "C2M":
            self.geo = "PLANE"
            self.qnd_mat = False
            self.reduce = self.c2m_reduce
            self.c2m = True

        if not is_bversion(410):
            self.shade_smooth = False

        # Mouse progress meter: Fake! just to show something is happening...
        # Actual (cheap+simple) progress tracking in console window std print out
        sys.stdout.write("\n[------------------ keI2M ---------------------]\n")
        if self.use_rgb:
            sys.stdout.write(" - Using RGB instead of Alpha -\n")
        self.wm = context.window_manager
        self.wm.progress_begin(0, 99)
        self.t = time.time()

        # Get current collection
        cvl = context.view_layer.active_layer_collection.name
        collections = [c.name for c in bpy.data.collections]
        if cvl in collections:
            self.coll = bpy.data.collections[cvl]
        else:
            self.coll = context.scene.collection

        # Batch setup ( since it refuses to reuse last used settings ran as op?..)
        if self.batch:
            self.opacity = k.opacity
            self.workres = k.workres
            self.geo = k.geo
            self.screw_flip = k.screw_flip
            self.screw_xcomp = k.screw_xcomp
            self.reduce = k.reduce
            self.shade_smooth = k.shade_smooth
            self.front_only = True
            self.qnd_mat = k.qnd_mat
            self.apply = k.apply
            self.apply_none = k.apply_none
            self.angle = k.angle
            self.vcolor = k.vcolor
            self.custom_workres = k.custom_workres

        # Auto Set View mode QoL (and make sure no geo smoothing is used for vertex color mode)
        if self.vcolor:
            self.reduce = "NONE"
            if context.space_data.shading.type in {"SOLID", "RENDERED", "MATERIAL"}:
                context.space_data.shading.color_type = 'VERTEX'
        else:
            if context.space_data.shading.type in {"SOLID", "RENDERED", "MATERIAL"}:
                context.space_data.shading.color_type = 'TEXTURE'

        # ----------------------------------------------------------------------------------------------
        # Image(s) Setup
        # ----------------------------------------------------------------------------------------------
        axis = ["Front", "Right", "Top", "Back", "Left", "Bottom"]
        front_img, right_img, top_img, back_img, left_img, bottom_img = None, None, None, None, None, None
        replace_right, replace_top, replace_left, replace_back, replace_bottom = False, False, False, False, False

        if k.FRONT:
            front_img = bpy.data.images.get(k.FRONT)
        if k.RIGHT:
            right_img = bpy.data.images.get(k.RIGHT)
        else:
            replace_right = True
        if k.TOP:
            top_img = bpy.data.images.get(k.TOP)
        else:
            replace_top = True
        if k.BACK:
            back_img = bpy.data.images.get(k.BACK)
        if k.LEFT:
            left_img = bpy.data.images.get(k.LEFT)
        else:
            replace_left = True
        if k.BOTTOM:
            bottom_img = bpy.data.images.get(k.BOTTOM)
        else:
            replace_bottom = True

        images = [front_img, right_img, top_img, back_img, left_img, bottom_img]

        count_images = [i for i in images if i is not None]
        if not count_images:
            self.report({"WARNING"}, "Aborted: Images not loaded in blend file")
            return {"CANCELLED"}

        # Missing Alpha Check
        if not alpha_check(count_images, self.use_rgb, self.c2m):
            self.report({"ERROR"}, "Aborted: Alpha channel missing!")
            return {"CANCELLED"}

        self.img_count = len(count_images)
        if self.img_count == 1 and self.geo == "BOOLEAN":
            self.geo = "PLANE"

        if self.img_count == 1 and self.geo != "BOOLEAN":
            self.front_only = True

        replaced_images = images.copy()

        if not self.front_only:
            if replace_right and replace_left:
                replaced_images[1] = front_img
                replaced_images[4] = front_img
            if replace_right and not replace_left:
                replaced_images[1] = left_img
            if replace_left and not replace_right:
                replaced_images[4] = right_img
            if replace_top and replace_bottom and not replace_right:
                replaced_images[2] = right_img
                replaced_images[5] = right_img
            elif replace_top and replace_bottom and replace_right:
                replaced_images[2] = front_img
                replaced_images[5] = front_img
            elif replace_bottom and not replace_top:
                replaced_images[2] = top_img
            elif replace_bottom:
                replaced_images[2] = front_img

        # ----------------------------------------------------------------------------------------------
        # Image resolution / square match check
        # ----------------------------------------------------------------------------------------------
        res_check = True
        ref_image = count_images[0]
        base_xy = ref_image.size
        if base_xy[1] != base_xy[0]:
            res_check = False
        else:
            bx = base_xy[0]
            for i in images:
                if i is not None:
                    x = i.size[0]
                    y = i.size[1]
                    if x != bx or y != bx or x != y:
                        res_check = False

        non_square_x = 1
        non_square_z = 1

        if not res_check:
            non_square_x = base_xy[0] / base_xy[1]
            non_square_z = base_xy[1] / base_xy[0]
            if non_square_z < non_square_x:
                non_square_x = 1
            else:
                non_square_z = 1

            if self.geo == "BOOLEAN":
                self.report({"INFO"}, "Warning: Images are not the same resolution -or- not square")

        # Work Res
        if self.workres == "IMAGE":
            work_res = ref_image.size[0]
            if work_res >= 1024:
                sys.stdout.write("WARNING: Work Resolution >= 1024 : May be slow/fail!\n")
        else:
            work_res = int(self.workres)

        if self.custom_workres != 0:
            work_res = self.custom_workres

        # ----------------------------------------------------------------------------------------------
        # Set Scale (& custom override)
        # ----------------------------------------------------------------------------------------------
        scl = 1 / work_res

        if self.size == "PIXEL":
            scl = self.pixel_width
            self.width = self.pixel_width * work_res
        elif self.size == "WIDTH":
            scl = self.width / work_res
            self.pixel_width = scl
        else:
            self.width = 1
            self.pixel_width = 0

        w = (work_res * scl) * 0.5

        # ----------------------------------------------------------------------------------------------
        # Process Mesh Images (0-2)
        # ----------------------------------------------------------------------------------------------

        # Main Object Name - w/o filetype
        name = ref_image.name
        if "." in name:
            name = name.split(".")[0]

        if self.geo == "BOOLEAN":
            mesh_images = images[:3]
            mesh_axis = axis[:3]
        else:
            images = [images[0]]
            axis = [axis[0]]
            mesh_images = images
            mesh_axis = axis

        if self.front_only and len(images) > 1:
            images = [front_img]
            axis = ["Front"]

        objects = []

        for image, axis_name in zip(mesh_images, mesh_axis):
            if image is not None:
                sys.stdout.write("%s Component:\n" % axis_name)
                img = image.copy()
                img.scale(width=work_res, height=work_res)

                # Make Pixel Map
                self.progress_update(context, " Generate Pixel Map ", False)
                pixels = img.pixels[:]
                width = img.size[0]
                height = img.size[1]
                pixel_map = self.make_pixel_map(width, height, pixels, scl, use_rgb=self.use_rgb)
                self.progress_update(context, " Generate Pixel Map ", True)

                # Create Mesh Data
                self.progress_update(context, " Create Mesh Data   ", False)
                mesh_name = name + "_i2m_" + axis_name
                existing = bpy.data.meshes.get(mesh_name)
                if existing:
                    bpy.data.meshes.remove(existing)
                mesh = self.make_mesh_data(pixel_map, work_res, scl, name=mesh_name, axis=axis_name)
                self.progress_update(context, " Create Mesh Data   ", True)

                # Bmesh Cleanup & Processing
                self.progress_update(context, " Mesh Cleanup       ", False)
                self.cleanup(mesh, scl, axis_name)

                # Create New Object from Mesh Data
                obj = self.make_scene_object(mesh, name=mesh_name)
                obj.data.uv_layers.new(name="UVmap")
                objects.append(obj)
                if axis_name == "Top":
                    obj.rotation_euler[2] = 1.5707963

                self.progress_update(context, " Mesh Cleanup       ", True)
                # Remove temp work image
                bpy.data.images.remove(img)

        final_object = objects.pop(0)

        sys.stdout.write("Object:\n")

        # ----------------------------------------------------------------------------------------------
        # # Make UV Projectors & Materials
        # ----------------------------------------------------------------------------------------------
        if self.c2m:
            for i, material in enumerate(self.cmats):
                mat = bpy.data.materials.new(name="I2M_"+str(i))
                mat.diffuse_color = material
                final_object.data.materials.append(mat)

        elif not self.vcolor:
            self.progress_update(context, " UV & Shading       ", False)

            projectors = []
            # Make UV Projectors
            if self.front_only:
                plist = ["Front"]
            else:
                plist = ["Front", "Right", "Top", "Back", "Left", "Bottom"]

            for axis_name in plist:
                pname = final_object.name + axis_name + "_UVProjector"
                projector = self.make_projector(pname, w, axis_name)
                projector.parent = final_object
                projectors.append(projector)

            mat_axis = []
            # Make Materials
            for img, axis_name in zip(images, axis):
                if img is not None:
                    # Adding Materials
                    material_name = img.name.split(".")[0] + "_Material"
                    mat = self.make_material(material_name, img)
                    final_object.data.materials.append(mat)
                    mat_axis.append(axis_name)

            if "Top" not in mat_axis and "Bottom" not in mat_axis:
                self.noz = True
            else:
                self.noz = False

        self.progress_update(context, " UV & Shading       ", True)

        # ----------------------------------------------------------------------------------------------
        # Modfiers & Finalize
        # ----------------------------------------------------------------------------------------------
        self.progress_update(context, " Modifiers/Finalize ", False)

        if self.geo == "SCREW":
            screw = final_object.modifiers.new(name="I2M Screw", type="SCREW")
            screw.use_merge_vertices = True

        elif self.geo == "BOOLEAN":
            solidify_objects = [final_object] + objects
            for obj in solidify_objects:
                context.view_layer.objects.active = obj
                solidify = obj.modifiers.new(name="I2M Solidify", type="SOLIDIFY")
                solidify.thickness = self.width
                solidify.offset = 0.0
                solidify.use_quality_normals = False
                if not self.apply_none:
                    bpy.ops.object.modifier_apply(modifier="I2M Solidify")

            final_object.select_set(True)
            context.view_layer.objects.active = final_object
            # BOOLEAN OPS
            for obj in objects:
                boolean = final_object.modifiers.new(name="I2MBool", type="BOOLEAN")
                boolean.operation = "INTERSECT"
                boolean.use_hole_tolerant = True
                boolean.object = obj
                if not self.apply_none:
                    bpy.ops.object.modifier_apply(modifier="I2MBool")
                    bpy.data.objects.remove(obj)

            # Assing materials to faces
            if not self.vcolor:
                idx, vecs = self.sort_material_slots(mat_axis)

                for p in final_object.data.polygons:
                    n = p.normal
                    p.material_index = 0
                    for i, v in zip(idx, vecs):
                        if -self.angle > n.dot(v) < self.angle:
                            p.material_index = i
                            break

        if self.shade_smooth:
            values = [True] * len(final_object.data.polygons)
            final_object.data.polygons.foreach_set("use_smooth", values)

        # Lastly, setup UV projection
        if self.vcolor or self.c2m:
            pass
        else:
            uv_project = final_object.modifiers.new(name="I2M UV-Project", type="UV_PROJECT")
            uv_project.uv_layer = "UVMap"
            uv_project.projector_count = len(projectors)
            for i, p in enumerate(projectors):
                uv_project.projectors[i].object = p
                if not self.front_only:
                    if i == 4 and self.flip_left:
                        p.scale = (-w, w, w)
                    elif i == 3 and self.flip_back:
                        p.scale = (-w, w, w)
            if self.noz and not self.front_only:
                uv_project.projectors[2].object = None
                uv_project.projectors[5].object = None

        # once more to make sure I guess
        final_object.select_set(state=True)
        context.view_layer.objects.active = final_object

        if self.apply:
            for m in final_object.modifiers:
                bpy.ops.object.modifier_apply(modifier=m.name)
            for p in projectors:
                bpy.data.objects.remove(p)

        if not res_check and not self.geo == "SCREW":
            final_object.scale.x = non_square_x
            final_object.scale.z = non_square_z

        self.progress_update(context, " Modifiers/Finalize ", True)

        self.wm.progress_end()

        if self.tot > 60:
            tot = str(round((self.tot / 60), 1)) + "min"
        else:
            tot = '{:f}'.format(self.tot).rstrip('0') + "s"

        sys.stdout.write("\n Total              : [   COMPLETE   ] %s\n\n" % tot)
        sys.stdout.flush()

        if not self.batch:
            k.opacity = self.opacity
            k.workres = self.workres
            k.geo = self.geo
            k.screw_flip = self.screw_flip
            k.screw_xcomp = self.screw_xcomp
            k.reduce = self.reduce
            k.shade_smooth = self.shade_smooth
            k.front_only = self.front_only
            k.qnd_mat = self.qnd_mat
            k.apply = self.apply
            k.apply_none = self.apply_none
            k.angle = self.angle
            k.custom_workres = self.custom_workres
            k.vcolor = self.vcolor

        # Needed for 1st-runs, or images can't be accessed by redo panel?!
        bpy.ops.ed.undo_push()
        context.area.tag_redraw()

        # Very silly
        if self.geo != "BOOLEAN" and self.front_only:
            self.front_only = False

        return {'FINISHED'}


class KeI2Mreload(Operator):
    bl_idname = "ke.i2m_reload"
    bl_label = "I2M Reload"
    bl_description = "Reload I2M image(s)"

    def execute(self, context):
        k = context.scene.kei2m
        front_img = bpy.data.images.get(k.FRONT)
        right_img = bpy.data.images.get(k.RIGHT)
        top_img = bpy.data.images.get(k.TOP)
        back_img = bpy.data.images.get(k.BACK)
        left_img = bpy.data.images.get(k.LEFT)
        bottom_img = bpy.data.images.get(k.BOTTOM)
        images = [front_img, right_img, top_img, back_img, left_img, bottom_img]
        kindex = ["FRONT", "RIGHT", "TOP", "BACK", "LEFT", "BOTTOM"]
        missing = 0
        found = 0
        for index, i in enumerate(images):
            if i is None:
                k[kindex[index]] = ""
                missing += 1
            if i is not None:
                i.reload()
        context.area.tag_redraw()
        if missing:
            self.report({"INFO"}, "%s slot(s) reloaded - %s slot(s) empty/invalid" % (str(found), str(missing)))
        return {"FINISHED"}


class KeI2Mfilebrowser(Operator, ImportHelper):
    bl_idname = "ke.i2m_filebrowser"
    bl_label = "Load Image(s)"
    bl_description = "Open filebrowser to load image(s)"

    filter_glob: StringProperty(
        default='*.png;*.tif;*.tiff;*.exr;*.hdr;*.tga;*.sgi;*.rgb;*.bw;*.jp2;*.j2c;*.cin;*.dpx',
        options={'HIDDEN'})

    axis : StringProperty(name="Axis Name", default="", options={'HIDDEN'})

    def execute(self, context):
        print("\n[------------- keI2M Image Loader -------------]")
        k = context.scene.kei2m
        kap = context.preferences.addons["ke_i2m"].preferences

        loaded = bpy.path.basename(self.filepath).split(".")
        loaded_name = loaded[0]

        # filename check
        if len(loaded_name) > 63:
            print("Filename: '%s'" % loaded_name)
            print("  -> Filename is too long (>64) - will be shortened (internally in Blender)\n"
                  "  -> 'Reload Image(s)' operator will not work\n")

        ext = "." + loaded[-1]
        loaded_dir = os.path.dirname(self.filepath)

        paths = [self.filepath]
        slots = [self.axis]

        if k.autofill and k.geo == "BOOLEAN":
            suffix = ["_front", "_right", "_top", "_back", "_left", "_bottom"]
            suffix_check = False

            if suffix[0] not in loaded_name:
                self.report({"INFO"}, "Failed to autoload %s - Wrong suffix for 1st slot" % loaded_name)
                return {"CANCELLED"}

            for s in suffix:
                if s in loaded_name:
                    suffix_check = True
                    loaded_name = loaded_name.replace(s, '')
                    suffix.remove(s)
                    break
            if suffix_check:
                for s in suffix:
                    n = loaded_name + s + ext
                    f = os.path.join(loaded_dir, n)
                    paths.append(f)
                    su = s.upper()[1:]
                    slots.append(su)
            else:
                print("Autoload failed - Invalid suffix")

        alpha_missing = []
        c2m = True if k.geo == "C2M" else False

        for f, slot in zip(paths, slots):
            img = load_slot(f)
            if img is not None:
                filename = img.name
                k[slot] = filename
                context.area.tag_redraw()
                if filename:
                    if alpha_check([img], kap.use_rgb, c2m):
                        print("Image Loaded: %s" % filename)
                    else:
                        alpha_missing.append(filename)
                else:
                    self.report({"WARNING"}, "File Not Loaded")
        if alpha_missing:
            self.report({"WARNING"}, "Missing Alpha Channel: %s" % alpha_missing)

        return {'FINISHED'}


class KeI2Mbatchbrowser(Operator, ImportHelper):
    bl_idname = "ke.i2m_batchbrowser"
    bl_label = "Batch Process Folder"
    bl_description = "Pick a folder to batch kei2m on ALL the images in the folder.\n --> Using last used settings <--"

    filter_glob : bpy.props.StringProperty(subtype="DIR_PATH")

    def execute(self, context):
        if not self.filepath:
            return {"CANCELLED"}

        # Load images in folder to Batch
        filter_glob = ('.png', '.tif', '.tiff', '.exr', '.hdr', '.tga', '.sgi', '.rgb', '.bw', '.jp2', '.j2c', '.cin',
                       '.dpx')
        images = []
        img_count = 0

        for file in os.listdir(self.filepath):
            if file.lower().endswith(filter_glob):
                path = os.path.join(self.filepath, file)
                img = load_slot(path)
                if img is not None:
                    img_count += 1
                    images.append(file)

        if not images:
            sys.stdout.write("\nkei2m Batch Process Aborted: No images could be loaded\n")
            self.report({"INFO"}, "Aborted: No images could be loaded")
            return {"CANCELLED"}

        sys.stdout.write("\nkei2m Batch Process Images Loaded: %s\n" % str(img_count))

        # Batch all loaded images
        k = context.scene.kei2m
        bpy.ops.ke.i2m_clearslot(axis="ALL")

        for img in images:
            k.FRONT = img
            bpy.ops.ke.i2m(batch=True)

        k.FRONT = ""
        return {'FINISHED'}


class KeI2Mclearslot(Operator):
    bl_idname = "ke.i2m_clearslot"
    bl_label = "Clear i2m File Slot"

    axis : StringProperty(name="Axis Name", default="")

    @classmethod
    def description(cls, context, properties):
        if properties.axis == "ALL":
            return "Clear All Image Slots"
        else:
            return "Clear Image Slot"

    def execute(self, context):
        k = context.scene.kei2m
        if self.axis == "ALL":
            for axis in ["FRONT", "RIGHT", "TOP", "BACK", "LEFT", "BOTTOM"]:
                k[axis] = ""
        else:
            k[self.axis] = ""
        context.area.tag_redraw()
        return {'FINISHED'}


def load_slot(path):
    try:
        img = bpy.data.images.load(path, check_existing=True)
    except (OSError, IOError, Exception):
        img = None
    return img


# ------------------------------------------------------------------------------------------------------------
# UI
# ------------------------------------------------------------------------------------------------------------
class VIEW3D_PT_i2m(Panel):
    bl_label = 'keI2M v%s' % str(kei2m_version)
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'kei2m'

    def draw_header_preset(self, context):
        layout = self.layout
        layout.emboss = 'NONE'
        row = layout.row(align=False)
        row.prop(context.scene.kei2m, 'info', text="", icon="QUESTION", icon_only=True, emboss=False)
        row.separator(factor=0.5)

    def draw(self, context):
        k = context.scene.kei2m
        has_front = bool(k.FRONT)
        has_2 = any([bool(k.RIGHT), bool(k.TOP)])

        layout = self.layout
        col = layout.column(align=True)
        row = col.row(align=True)
        row.enabled = False
        row.alignment = "LEFT"
        row.label(text="Select Mode")
        row = col.row(align=True)
        row.prop(k, "geo", text="", expand=False)
        row = col.row(align=True)
        row.enabled = False
        row.alignment = "LEFT"

        if k.geo == "BOOLEAN":
            row.label(text="Load Images")
            row = layout.row(align=True)
            row.prop(k, "autofill")
            row = layout.row(align=True)
            row.operator("ke.i2m_filebrowser", text="Front   ", icon="FILEBROWSER").axis = "FRONT"
            row.prop(k, "FRONT", text="")
            row.operator("ke.i2m_clearslot", text="", icon="X").axis = "FRONT"
            row = layout.row(align=True)
            if not has_front:
                row.enabled = False
            row.operator("ke.i2m_filebrowser", text="Right   ", icon="FILEBROWSER").axis = "RIGHT"
            row.prop(k, "RIGHT", text="")
            row.operator("ke.i2m_clearslot", text="", icon="X").axis = "RIGHT"
            row = layout.row(align=True)
            if not has_front:
                row.enabled = False
            row.operator("ke.i2m_filebrowser", text="Top      ", icon="FILEBROWSER").axis = "TOP"
            row.prop(k, "TOP", text="")
            row.operator("ke.i2m_clearslot", text="", icon="X").axis = "TOP"

            row = layout.row(align=True)
            row.scale_y = 0.5
            row.enabled = False
            row.label(text="Extra Textures Only:")

            row = layout.row(align=True)
            if not has_front or not has_2:
                row.enabled = False
            row.scale_y = 1
            row.operator("ke.i2m_filebrowser", text="Back    ", icon="FILEBROWSER").axis = "BACK"
            row.prop(k, "BACK", text="")
            row.operator("ke.i2m_clearslot", text="", icon="X").axis = "BACK"
            row = layout.row(align=True)
            if not has_front or not has_2:
                row.enabled = False
            row.operator("ke.i2m_filebrowser", text="Left     ", icon="FILEBROWSER").axis = "LEFT"
            row.prop(k, "LEFT", text="")
            row.operator("ke.i2m_clearslot", text="", icon="X").axis = "LEFT"
            row = layout.row(align=True)
            if not has_front or not has_2:
                row.enabled = False
            row.operator("ke.i2m_filebrowser", text="Bottom", icon="FILEBROWSER").axis = "BOTTOM"
            row.prop(k, "BOTTOM", text="")
            row.operator("ke.i2m_clearslot", text="", icon="X").axis = "BOTTOM"
            row = layout.row(align=True)
            if not has_front:
                row.enabled = False
            row.operator("ke.i2m_clearslot", text="Clear All", icon="CANCEL").axis = "ALL"
        else:
            row.label(text="Load Image")
            row = col.row(align=True)
            row.operator("ke.i2m_filebrowser", text="", icon="FILEBROWSER").axis = "FRONT"
            row.prop(k, "FRONT", text="")
            row.operator("ke.i2m_clearslot", text="", icon="X").axis = "FRONT"

        row = layout.row(align=False)
        row.operator("ke.i2m_reload", text="Reload Image(s)", icon="LOOP_BACK")
        box = layout.box()
        box.operator("ke.i2m", text="Reset To Defaults").reset = True
        box.operator("ke.i2m_batchbrowser", icon="FILE_FOLDER")

        row = box.row()
        row.scale_y = 1.8
        if not has_front:
            row.enabled = False
        if k.geo == "BOOLEAN" and not has_2:
            row.enabled = False
        if k.geo == "BOOLEAN":
            row.operator("ke.i2m", text="Image To Mesh").reduce = "DISSOLVE"
        else:
            row.operator("ke.i2m", text="Image To Mesh")


# ------------------------------------------------------------------------------------------------------------
# Prefs
# ------------------------------------------------------------------------------------------------------------
# Panels to update
panels = (
        VIEW3D_PT_i2m,
        )


def update_panel(self, context):
    message = "kei2m : panel update failed"
    try:
        for panel in panels:
            if "bl_rna" in panel.__dict__:
                bpy.utils.unregister_class(panel)

        for panel in panels:
            panel.bl_category = context.preferences.addons[__name__].preferences.category
            bpy.utils.register_class(panel)

    except Exception as e:
        print("\n[{}]\n{}\n\nError:\n{}".format(__name__, message, e))
        pass


class KeI2Maddonprefs(AddonPreferences):
    bl_idname = __name__

    category: StringProperty(
            name="Tab Category",
            description="Choose a name (category) for tab placement",
            default="kei2m",
            update=update_panel
            )
    use_rgb: BoolProperty(name="Use RGB instead of Alpha", default=False,
                          description="Use the user RGB color instead of Alpha channel")
    user_rgb: FloatVectorProperty(name="User RGB", subtype="COLOR", size=3, default=(1.0, 1.0, 1.0))
    cap: IntProperty(default=16, name="Material Cap",
                     description="Maximum number of materials generated in Color 2 Material Mode.")

    def draw(self, context):
        layout = self.layout
        row = layout.row()
        row.label(text="Tab Location (Category):")
        row.prop(self, "category", text="")
        row = layout.row(align=True)
        row.prop(self, "use_rgb", toggle=True)
        row.prop(self, "user_rgb", text="")
        row = layout.row()
        row.use_property_split = True
        row.prop(self, "cap")


class KeI2Mprops(PropertyGroup):
    FRONT: StringProperty(default="", name="Front Image", description="-Y Axis Image")
    RIGHT: StringProperty(default="", name="Right Image", description="+X Axis Image")
    TOP: StringProperty(default="", name="Top Image", description="+Z Axis Image")
    BACK: StringProperty(default="", name="Back Image", description="+Y Axis Image")
    LEFT: StringProperty(default="", name="Left Image", description="-X Axis Image")
    BOTTOM: StringProperty(default="", name="Bottom Image", description="-Z Axis Image")
    autofill: BoolProperty(name="Autofill", default=False,
                           description="Auto-Fill other slots by suffixes when loading an image.\n "
                                       "Valid suffixes:\n"
                                       "_front, _right, _top\n"
                                       "_back,  _left,  _bottom")
    info: BoolProperty(name="kei2m General Info", default=False,
                       description="Load Alpha Images to convert into mesh (in Object Mode).\n"
                                   "A single image is required for PLANE & SCREW (& C2M) modes.\n"
                                   "2 or more images are required for BOOLEAN mode.\n"
                                   "Adjust result in Redo Panel.\n"
                                   "Non-alpha color option in Addon Preferences.\n"
                                   "See progress output in Console Window.")
    opacity: IntProperty(default=95)
    workres: StringProperty(default="128")
    geo: EnumProperty(
        items=[("PLANE", "Plane", "Converts into a Plane Mesh", "", 1),
               ("SCREW", "Screw", "A Screw Modifier revolves one side of a (symmetrical) image to a cylindrical shape",
                "", 2),
               ("BOOLEAN", "Boolean",
                "Uses 2-3 images (front + right/top) to carve out a rough 3d shape using Intersect", "", 3),
               ("C2M", "Color 2 Material", "Splits along color borders and assigns color materials on an mesh plane",
                "", 4),
               ],
        name="Geo Type", default="PLANE",
        description="Which type of geo conversion to use")
    screw_flip: BoolProperty(default=False)
    screw_xcomp: IntProperty(default=15)
    reduce: StringProperty(default="SIMPLE")
    shade_smooth: BoolProperty(default=True)
    front_only: BoolProperty(default=False)
    qnd_mat: BoolProperty(default=False)
    apply: BoolProperty(default=False)
    apply_none: BoolProperty(default=False)
    angle: FloatProperty(default=0.5)
    custom_workres: IntProperty(default=0)
    vcolor: BoolProperty(default=False)

# ------------------------------------------------------------------------------------------------------------
# Registration
# ------------------------------------------------------------------------------------------------------------
classes = (
    KeI2M,
    VIEW3D_PT_i2m,
    KeI2Maddonprefs,
    KeI2Mprops,
    KeI2Mfilebrowser,
    KeI2Mreload,
    KeI2Mclearslot,
    KeI2Mbatchbrowser,
)


def register():
    for c in classes:
        bpy.utils.register_class(c)

    bpy.types.Scene.kei2m = PointerProperty(type=KeI2Mprops)

    # Force Panel Udpate, for custom tab location.
    try:
        if "bl_rna" in VIEW3D_PT_i2m.__dict__:
            bpy.utils.unregister_class(VIEW3D_PT_i2m)
        VIEW3D_PT_i2m.bl_category = bpy.context.preferences.addons[__name__].preferences.category
        bpy.utils.register_class(VIEW3D_PT_i2m)
    except Exception as e:
        print("kei2m panel update failed:\n", e)


def unregister():
    for c in reversed(classes):
        bpy.utils.unregister_class(c)

    try:
        del bpy.types.Scene.kei2m

    except Exception as e:
        print('unregister fail:\n', e)
        pass


if __name__ == "__main__":
    register()
