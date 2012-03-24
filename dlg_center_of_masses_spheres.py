# -*- coding: utf-8 -*-
#
#   Copyright (C) 2011, 2012 Jan-Philip Gehrcke
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
#


from pymol import cmd
from collections import defaultdict
import itertools
import math
import operator
import os
import time
import numpy

def dlg_models_pdb_format(dlg_raw_lines):
    """Return all lines belonging to a model (in PDB(QT) format)
    """
    startmarker = "DOCKED: MODEL" # should work for AD3 and 4
    endmarker = "DOCKED: ENDMDL"
    modelblock = False
    modellines = []
    for line in dlg_raw_lines:
        if not modelblock:
            if line.startswith(startmarker):
                modelblock = True
        else:
            if line.startswith(endmarker):
                modelblock = False
                #print "len model: %s" % len(modellines)
                yield modellines
                modellines = []
        if modelblock:
            if len(line) > 8:
                # do not use the "DOCKED: " at the beginning of the line
                modellines.append(line[8:])


def dlg_models_to_objects(dlg_prefix, directory):
    """
    In current working directory, read all files starting with `dlg_prefix`.
    Assume that they are in Autodock (4) DLG format and extract ligand models
    as PDB strings. Create PyMOL objects from these PDB strings. Yield these
    objects.
    """
    if not directory:
        directory = os.getcwd()
    if not os.path.isdir(directory):
        print "not a directory: `%s`" % directory
        return
    for filename in os.listdir(directory):
        if filename.startswith(dlg_prefix):
            filepath = os.path.join(directory, filename)
            if not os.path.isfile(filepath):
                continue
            print "reading file: %s" % filename
            with open(filepath) as f:
                dlg_raw_lines = f.readlines()
            for modellines in dlg_models_pdb_format(dlg_raw_lines):
                p = os.path.splitext(os.path.basename(filename))[0]
                objname = cmd.get_unused_name(prefix=p)
                #print objname
                pdb_string = '\n'.join(modellines)
                cmd.read_pdbstr(pdb_string, objname)
                yield objname


def load_all_dlg_models(dlg_prefix, directory=None):
    """
    For all Autodock ligand models in all files starting with `dlg_prefix`
    create an object.
    """
    # read file(s), create objects, get all object names
    objects = [o for o in dlg_models_to_objects(dlg_prefix, directory)]
    cmd.center()
    cmd.zoom()


def coc_spheres_for_dlg_prefix(dlg_prefix, directory=None, radius=0.5):
    """
DESCRIPTION

    For all Autodock ligand models in all files starting with
    `dlg_prefix` draw a sphere at the center of coordinates.

USAGE

   coc_spheres_for_dlg_prefix dlg_prefix, [directory=CWD, [radius=0.5]]


EXAMPLES

   coc_spheres_for_dlg_prefix 0.5, lga, /data/docking

AUTHOR

    Jan-Philip Gehrcke
    TU Dresden
    """
    # read file(s), create objects, get all object names
    objects = [o for o in dlg_models_to_objects(dlg_prefix, directory)]

    if not objects:
        print "No models extracted from file(s)"
        return

    # create COC-spheres for all objects

    # put them all to one object in different states, then show all states
    n = cmd.get_unused_name(prefix="%s_COC" % dlg_prefix)
    for i, objname in enumerate(objects):
        cmd.pseudoatom(n, pos=coc(objname), state=i)
        cmd.delete(objname)
    cmd.show("spheres", n)
    cmd.set("all_states", True)
    # OpenGL shaders: http://www.pymolwiki.org/index.php/Spheres
    cmd.set("sphere_mode", 5)
    cmd.show_as("spheres", n)

    cmd.alter(n, "vdw=%s" % radius)
    cmd.set(name="sphere_color", value="green", selection=n)

    # the alternative way: each sphere is its own object
    #[mutate_object_to_coc_sphere(o) for o in objects]
    cmd.center()
    cmd.zoom()


def mutate_object_to_coc_sphere(objname):
    """
    Calculate center of coordinates for object with name `objname`. Draw
    sphere at this point (pseudoatom)
    """
    #p = "%s_COC" % objname
    #n= cmd.get_unused_name(prefix=p)
    #cmd.pseudoatom(n, pos=coc(objname))
    #cmd.delete(objname)

    n = "%s_COC" % objname
    cmd.pseudoatom(n, pos=coc(objname))
    cmd.delete(objname)
    cmd.show("spheres", n)


def coc(objname):
    """Return center of coordinates for object
    """
    # This is the way I found to access the atoms in a PyMOL "model"
    #print "objectname: %s" % objname
    m = cmd.get_model(objname)
    atomcount = len(m.atom)
    #print "atomcount: %s" % atomcount

    center_x = sum((a.coord[0] for a in m.atom))/atomcount
    center_y = sum((a.coord[1] for a in m.atom))/atomcount
    center_z = sum((a.coord[2] for a in m.atom))/atomcount

    return (center_x, center_y, center_z)

cmd.extend("coc_spheres_for_dlg_prefix", coc_spheres_for_dlg_prefix)
cmd.extend("load_all_dlg_models", load_all_dlg_models)


def filtered_sphere_histogram(centers, sphere_radius):
    """
    `centers`: tuple of coordinate lists (used as sphere centers)
    `sphere_rad`: sphere radius
    For each sphere center, the number of other sphere centers (ligand COCS)
    within the sphere will be determined and saved to a histogram `h`. Return
    most occupied, not overlapping sphere centers including occupancy.
    """
    #print "len(centers): %s" % len(centers)
    t = time.time()
    h = defaultdict(int)
    for center_i in centers:
        #print "i: %s ; center_i: %s" % (i, center_i)
        for center_j in centers:
            d = distance(center_i, center_j)
            #print "distance: %s, sphere_radius: %s" % (d, sphere_radius)
            #print "type(d): %s, type(radius): %s" % (type(d), type(sphere_radius))
            if d < sphere_radius:
                #print "ADD TO CENTER %s" % (center_i, )
                h[center_i] += 1
    print "h creation time (s): %s" % (time.time()-t, )
    t = time.time()
    sorted_histogram = sorted(
        h.iteritems(), key=operator.itemgetter(1), reverse=True)
    print "h sortion time (s): %s" % (time.time()-t, )
    #print "sorted:"
    #for center, occupancy in sorted_histogram:
    #    print "%s: %s" % (center, occupancy)
    t = time.time()
    filtered_histogram = uniquify_centers(sorted_histogram, sphere_radius)
    print "h filter time (s): %s" % (time.time()-t, )
    return filtered_histogram


def uniquify_centers(sorted_histogram, sphere_radius):
    """
    Only keep center coordinates defining a sphere that is not overlapping with
    another sphere. It is impor-
    tant to start with the sorted histogram, because then the sphere with the
    highest occupancies are kept.
    """
    kept_centers = []
    for center, occupancy in sorted_histogram:
        keep = True
        for kept_center in (c for c, o in kept_centers):
            if distance(center, kept_center) <= 2 * sphere_radius:
                keep = False
                break
        if keep:
            kept_centers.append((center, occupancy))
    return kept_centers


def ligand_centers(dlg_prefix, directory):
    # read file(s), create objects, get all object names
    t = time.time()
    objects = [o for o in dlg_models_to_objects(dlg_prefix, directory)]
    print "Object creation time (s): %s" % (time.time()-t, )
    if not objects:
        print "No models extracted from file(s)"
        return
    # get coordinate centers of ligands
    t = time.time()
    centers = [coc(o) for o in objects]
    print "COC calc time (s): %s" % (time.time()-t, )
    # delete ligands
    t = time.time()
    [cmd.delete(o) for o in objects]
    print "Object deletetion time (s): %s" % (time.time()-t, )
    return centers


def show_most_occupied_spherical_volumes_cmd(
        dlg_prefix, radius, first_n, directory=None, coordinate_file=None):
    """
    Print filtered sphere histogram, draw most occupied spheres.
    """
    centers = ligand_centers(dlg_prefix, directory)
    filtered_histogram = filtered_sphere_histogram(
        centers=centers, sphere_radius=float(radius))
    # OpenGL shaders: http://www.pymolwiki.org/index.php/Spheres
    cmd.set("sphere_mode", 5)
    highest_occupancy = filtered_histogram[0][1]

    # write coc coordinates to file so that each line contains three numbers
    # that can be directly inserted in the autogrid parameter file
    if coordinate_file:
        f = open(coordinate_file, "w")

    for center, occupancy in filtered_histogram[:int(first_n)]:
        n = cmd.get_unused_name(prefix="most_occupied_spheres")
        print "%s: %f,%f,%f: %s" % (n, center[0], center[1], center[2], occupancy)
        if coordinate_file:
            f.write("%s\n" % " ".join(["%f" % c for c in center]))
        cmd.pseudoatom(n, pos=center)
        cmd.show("spheres", n)
        r = float(radius) * (float(occupancy)/highest_occupancy)
        cmd.alter(n, "vdw=%s" % r)
        cmd.set(name="sphere_color", value="red", selection=n)
    cmd.center()
    cmd.zoom()

    if coordinate_file:
        f.close()


def distance(x, y):
    return math.sqrt(sum(((e1-e2)**2 for e1, e2 in zip(x, y))))


cmd.extend("show_most_occupied_spherical_volumes", show_most_occupied_spherical_volumes_cmd)

