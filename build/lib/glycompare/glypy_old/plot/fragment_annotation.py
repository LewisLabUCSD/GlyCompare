from matplotlib import patches
from matplotlib import transforms as mtransforms
from matplotlib.textpath import TextPath
from glypy.plot.symbolic_nomenclature import line_to


def isclose(a, b):
    return abs(a - b) < 1e-4


rot90 = mtransforms.Affine2D().rotate_deg(90)


class BondCleavageArtist(object):

    def __init__(self, draw_tree, fragment, ax, label=True, **graphic_options):
        self.draw_tree = draw_tree
        self.fragment = fragment
        self.ax = ax
        self.artists = []
        self.patches = []
        self.graphic_options = graphic_options
        self.handle_glycosidic()
        self.handle_crossring()
        self.propagate_transform()

    def handle_glycosidic(self):
        for link_id, frag_type in self.fragment.link_ids.items():
            parent, child = self.draw_tree.get_link_pair(link_id)
            if frag_type[1] in "YZ":
                reducing = True
            else:
                reducing = False
            gba = GlycosidicBondCleavageStroke(self.draw_tree, parent, child, reducing,
                                               ax=self.ax, fragment=self.fragment,
                                               **self.graphic_options)
            gba.draw()
            self.artists.append(gba)
            self.patches.extend(gba.patches)

    def handle_crossring(self):
        for crossring_id, frag_type in self.fragment.crossring_cleavages.items():
            node = self.draw_tree.get(crossring_id)
            if frag_type[1] == "X":
                reducing = True
            else:
                reducing = False
            cra = CrossRingBondCleavageStroke(
                self.draw_tree, node, reducing, self.ax, self.fragment,
                **self.graphic_options)
            cra.draw()
            self.artists.append(cra)
            self.patches.extend(cra.patches)

    def propagate_transform(self):
        transform = self.draw_tree.transform
        if transform is None:
            return
        transform = transform + self.ax.transData
        for patch in self.patches:
            patch.set_transform(transform)


class CleavageStrokeArtistBase(object):

    def __init__(self, draw_tree, fragment, reducing, ax, **options):
        self.draw_tree = draw_tree
        self.ax = ax
        self.reducing = reducing
        self.fragment = fragment
        self.patches = []
        self.x_coords = []
        self.y_coords = []
        self.patch_dict = draw_tree.data['patches']
        self.options = options
        self.edge_size = None

    def save_patch(self, key, value):
        self.patches.append(value)
        if key in self.patch_dict:
            self.patch_dict[key].append([value])
        else:
            self.patch_dict[key] = [[value]]

    def textpath(self, x, y, text, line_weight=0.5, **kwargs):
        fs = kwargs.get("fontsize", 2) * .5
        t_path = TextPath((x, y), s=text, size=fs)
        patch = patches.PathPatch(t_path, facecolor="black", lw=line_weight / 20., zorder=4)
        a = self.ax.add_patch(patch)
        return a


class CrossRingBondCleavageStroke(CleavageStrokeArtistBase):

    def __init__(self, draw_tree, node, reducing, ax, fragment, **options):
        super(CrossRingBondCleavageStroke, self).__init__(
            draw_tree, fragment, reducing, ax, **options)
        self.node = node
        self.mode = None

    def draw(self):
        self.compute_positions()
        self.draw_main_stroke()
        self.draw_edge_stroke()

    def compute_positions(self):
        px, py = self.node.x, self.node.y
        if self.node.children:
            nx, ny = -1, -1
            for next_node in self.node.children:
                cx, cy = next_node.x, next_node.y
                if cy > ny:
                    nx, ny = cx, cy
            self.mode = 'child'
            cx, cy = nx, ny
            xcenter = (max(px, cx) - min(px, cx)) / 2. + min(px, cx)
            ycenter = (max(py, cy) - min(py, cy)) / 2. + min(py, cy)
            xlength = (xcenter - min(px, cx)) / 2.
            ylength = (ycenter - min(py, cy)) / 2.
            if py < cy:
                lx2, ly2 = self.node.x + xlength, self.node.y + ylength
                lx1, ly1 = self.node.x - xlength, self.node.y - ylength
            else:
                lx1, ly1 = self.node.x - xlength, self.node.y + ylength
                lx2, ly2 = self.node.x + xlength, self.node.y - ylength
            edge_size = abs(ly1 - ly2)

        elif self.node.parent:
            nx, ny = self.parent
            self.mode = 'parent'

            cx, cy = px, py
            px, py = nx, ny

            xcenter = (max(px, cx) - min(px, cx)) / 2. + min(px, cx)
            ycenter = (max(py, cy) - min(py, cy)) / 2. + min(py, cy)
            xlength = (xcenter - min(px, cx)) / 2.
            ylength = (ycenter - min(py, cy)) / 2.
            if py < cy:
                lx2, ly2 = xcenter + xlength, ycenter + ylength
                lx1, ly1 = xcenter - xlength, ycenter - ylength
            else:
                lx1, ly1 = xcenter - xlength, ycenter + ylength
                lx2, ly2 = xcenter + xlength, ycenter - ylength
            edge_size = abs(ly1 - ly2)

        else:
            print("Orphaned Node. Cross-ring fragment drawing may be bad")

        self.x_coords = [lx1, lx2]
        self.y_coords = [ly1, ly2]
        self.edge_size = edge_size

    def draw_main_stroke(self, **kwargs):
        graphic_options = self.options.copy()
        graphic_options.update(kwargs)

        color = graphic_options.get('color', 'red')
        lw = graphic_options.get('lw', 2)

        lx1, lx2 = self.x_coords
        ly1, ly2 = self.y_coords

        line = line_to(self.ax, lx1, ly1, lx2, ly2,
                       color=color, zorder=10, lw=lw)
        self.save_patch(self.fragment.name, line)
        line.set_gid(self.draw_tree.uuid + '-' + self.fragment.name)

    def draw_edge_stroke(self, **kwargs):
        graphic_options = self.options.copy()
        graphic_options.update(kwargs)

        color = graphic_options.get('color', 'red')
        lw = graphic_options.get('lw', 2)

        lx1, lx2 = self.x_coords
        ly1, ly2 = self.y_coords

        if self.reducing:
            line = line_to(self.ax, lx1, ly2, lx1, ly2 -
                           self.edge_size, color=color, zorder=3, lw=lw)
            line.set_gid(self.draw_tree.uuid + '-' +
                         self.fragment.name + "-direction")
            self.save_patch(self.fragment.name + "_direction", line)
        else:
            line = line_to(self.ax, lx2, ly2, lx2, ly2 +
                           self.edge_size, color=color, zorder=3, lw=lw)
            line.set_gid(self.draw_tree.uuid + '-' +
                         self.fragment.name + "-direction")
            self.save_patch(self.fragment.name + "_direction", line)


class GlycosidicBondCleavageStroke(CleavageStrokeArtistBase):

    def __init__(self, draw_tree, parent, child, reducing, ax, fragment, **options):
        super(GlycosidicBondCleavageStroke, self).__init__(
            draw_tree, fragment, reducing, ax, **options)

        self.parent = parent
        self.child = child

        self.mode = None

    def draw(self):
        self.compute_positions()
        self.draw_main_stroke()
        self.draw_edge_stroke()

    def compute_positions(self):
        parent, child = self.parent, self.child
        px, py = parent.x, parent.y
        cx, cy = child.x, child.y

        if isclose(py, cy) and not isclose(px, cx):
            mode = 'offsequence'
            center = (max(px, cx) - min(px, cx)) / 2. + min(px, cx)
            length = (center - min(px, cx))
            lx1, ly1 = center, py + length
            lx2, ly2 = center, py - length
            edge_size = abs(ly1 - ly2) / 2.

        elif not isclose(py, cy) and isclose(px, cx):
            mode = "sequential"
            center = (max(py, cy) - min(py, cy)) / 2. + min(py, cy)
            length = center - min(py, cy)
            lx1, ly1 = px + length, center
            lx2, ly2 = px - length, center
            edge_size = abs(lx1 - lx2) / 2.
        else:
            mode = "branching"
            xcenter = (max(px, cx) - min(px, cx)) / 2. + min(px, cx)
            ycenter = (max(py, cy) - min(py, cy)) / 2. + min(py, cy)
            xlength = (xcenter - min(px, cx)) / 2.
            ylength = (ycenter - min(py, cy)) / 2.
            if py < cy:
                lx2, ly2 = xcenter + xlength, ycenter + ylength
                lx1, ly1 = xcenter - xlength, ycenter - ylength
            else:
                lx1, ly1 = xcenter - xlength, ycenter + ylength
                lx2, ly2 = xcenter + xlength, ycenter - ylength
            edge_size = abs(ly1 - ly2)

        self.mode = mode
        self.edge_size = edge_size
        self.x_coords = [lx1, lx2]
        self.y_coords = [ly1, ly2]

    def draw_main_stroke(self, **kwargs):
        graphic_options = self.options.copy()
        graphic_options.update(kwargs)

        color = graphic_options.get('color', 'red')
        lw = graphic_options.get('lw', 2)

        lx1, lx2 = self.x_coords
        ly1, ly2 = self.y_coords

        line = line_to(self.ax, lx1, ly1, lx2, ly2,
                       color=color, zorder=3, lw=lw)
        self.save_patch(self.fragment.name, line)
        line.set_gid(self.draw_tree.uuid + '-' + self.fragment.name)

    def draw_edge_stroke(self, **kwargs):
        graphic_options = self.options.copy()
        graphic_options.update(kwargs)

        color = graphic_options.get('color', 'red')
        lw = graphic_options.get('lw', 2)

        lx1, lx2 = self.x_coords
        ly1, ly2 = self.y_coords
        # if self.mode == "offsequence":
        #    lx1, lx2, _ = rot90.get_matrix().dot([lx1, lx2, 0])
        #    ly1, ly2, _ = rot90.get_matrix().dot([ly1, ly2, 0])
        if self.reducing:
            line = line_to(self.ax, lx1, ly2, lx1, ly2 -
                           self.edge_size, color=color, zorder=3, lw=lw)
            line.set_gid(self.draw_tree.uuid + '-' +
                         self.fragment.name + "-direction")
            self.save_patch(self.fragment.name + "_direction", line)
        else:
            line = line_to(self.ax, lx2, ly2, lx2, ly2 +
                           self.edge_size, color=color, zorder=3, lw=lw)
            line.set_gid(self.draw_tree.uuid + '-' +
                         self.fragment.name + "-direction")
            self.save_patch(self.fragment.name + "_direction", line)


# Old Reference Implementation
def draw_cleavage(self, fragment=None, at=(0, 0), ax=None, scale=0.1, color='red', label=True):  # pragma: no cover
    '''
    .. warning::
        Here be magical numbers. This method is highly experimental and may not behave nicely
        when projected onto a transformed topology

    '''
    if ax is None:
        ax = self.axes
    if ax is None:
        raise ValueError("`ax` is required")
    if fragment is None:
        raise ValueError("`fragment` is required")

    break_targets = fragment.link_ids
    crossring_targets = fragment.crossring_cleavages

    def isclose(a, b):
        return abs(a - b) < 1e-4

    def save_patch(data, key, value):
        _created_patches.append(value)
        if key in data:
            data[key].append([value])
        else:
            data[key] = [[value]]

    def get_scale(trans_mat):
        ty = trans_mat[1, 0]
        return ty

    def textpath(x, y, text, line_weight=0.5, **kwargs):
        fs = kwargs.get("fontsize", 2) * scale * .5
        t_path = TextPath((x, y), s=text, size=fs)
        patch = patches.PathPatch(t_path, facecolor="black", lw=line_weight / 20., zorder=4)
        a = ax.add_patch(patch)
        return a

    if self.transform:
        trans_scale = max(get_scale(self.transform.get_matrix()) * 3e-5, 0.5)
    else:
        trans_scale = 1

    scale *= 2 * trans_scale
    print(scale, trans_scale)

    patch_dict = self.data['patches']
    _created_patches = []

    for link_break in break_targets:
        parent, child = self.get_link_pair(link_break)
        px, py = parent.coords(at)
        cx, cy = child.coords(at)

        print('p', (px, py), 'c', (cx, cy))
        branch_point = False
        length = scale
        if isclose(py, cy) and not isclose(px, cx):
            center = (max(px, cx) - min(px, cx)) / 2. + min(px, cx)
            length = max((center - min(px, cx)) * trans_scale, scale)
            if length < 1:
                length = length / length * 0.25
            else:
                length = length * length * 0.5
            lx1, ly1 = center, py + length
            lx2, ly2 = center, py - length
            edge_size = max(abs(ly1 - ly2) / 4000., scale * 2)
        elif not isclose(py, cy) and isclose(px, cx):
            center = (max(py, cy) - min(py, cy)) / 2. + min(py, cy)
            length = max((center - min(py, cy)) * trans_scale, scale)
            if length < 1:
                length = length / length * 0.25
            else:
                length = length * length * 0.5
            lx1, ly1 = px + length, center
            lx2, ly2 = px - length, center
            edge_size = max(abs(lx1 - lx2) / 4000., scale * 2)
        else:
            branch_point = True
            xcenter = (max(px, cx) - min(px, cx)) / 2. + min(px, cx)
            ycenter = (max(py, cy) - min(py, cy)) / 2. + min(py, cy)
            xlength = (xcenter - min(px, cx)) * trans_scale
            ylength = (ycenter - min(py, cy)) * trans_scale
            if py < cy:
                lx2, ly2 = xcenter + scale, ycenter + scale
                lx1, ly1 = xcenter - scale, ycenter - scale
            else:
                lx1, ly1 = xcenter - scale, ycenter + scale
                lx2, ly2 = xcenter + scale, ycenter - scale
            edge_size = scale * 2  # max(min(abs(lx1 - lx2) / 40000., abs(ly1 - ly2) / 40000.),scale * 2)

        print("Main Line: ", (lx1, lx2), (ly1, ly2, ly2 - ly1))

        lw = 2

        line = ax.plot((lx1, lx2), (ly1, ly2), color=color, zorder=3, lw=lw)
        save_patch(patch_dict, fragment.name, line[0])
        self.data['position'][fragment.name] = (lx1, lx2), (ly1, ly2)
        line[0].set_gid(self.uuid + '-' + fragment.name)

        print((lx1, lx2), (ly1, ly2), length, edge_size)

        if fragment.link_ids[link_break][1] in {"B", "C"}:
            label_x = (lx2, lx2 - edge_size)
            if branch_point:
                label_x = [x - 2 * scale for x in label_x]

            line = ax.plot(label_x, (ly1, ly1), color=color, zorder=3, lw=lw)
            line[0].set_gid(self.uuid + '-' + fragment.name + "-direction")
            save_patch(patch_dict, fragment.name + "_direction", line[0])
            self.data['position'][fragment.name + "_direction"] = label_x, (ly1, ly1)
            if label:
                text = textpath(label_x[0] - .4, ly1 + 0.05, fragment.fname)
                save_patch(patch_dict, fragment.name + "_text", text)
                self.data['position'][fragment.name + "_text"] = label_x[0] - .4, ly1 + 0.05
        else:
            line = ax.plot((lx2, lx2 + edge_size), (ly2, ly2), color=color, zorder=3, lw=lw)
            line[0].set_gid(self.uuid + '-' + fragment.name + "-direction")
            save_patch(patch_dict, fragment.name + "_direction", line[0])
            self.data['position'][fragment.name + "_direction"] = (lx2, lx2 + scale), (ly2, ly2)

            if label:
                lx2 += (0.05 * scale)
                ly2 -= (0.15 * scale)
                save_patch(patch_dict, fragment.name + "_text", textpath(lx2 + 0.05, ly2 - 0.15,
                                                                        fragment.fname))
                self.data['position'][fragment.name + "_text"] = lx2 + 0.05, ly2 - 0.15

    for crossring in crossring_targets:
        target = self.get(crossring)
        cx, cy = target.coords(at)
        line = ax.plot((cx - scale, cx + scale), (cy + scale, cy - scale), color=color, zorder=3)
        patch_dict[fragment.name] = line[0]
        self.data['position'][fragment.name] = (cx - scale, cx + scale), (cy + scale, cy - scale)
        line[0].set_gid(self.uuid + '-' + fragment.name.replace(",", '_'))
        annotation_name = re.sub(r'\d,\d', '', fragment.fname)
        if fragment.crossring_cleavages[crossring][1] == "X":
            line = ax.plot((cx + scale, cx + 2 * scale), (cy - scale, cy - scale), color=color, zorder=3)
            line[0].set_gid(self.uuid + '-' + fragment.name + "-direction")
            patch_dict[fragment.name + "_direction"] = line[0]
            self.data['position'][fragment.name + "_direction"] =\
                (cx + scale, cx + 2 * scale), (cy - scale, cy - scale)

            if label:
                ax.text((cx + scale) - 0.4, (cy - scale) - .15, annotation_name)
        else:
            line = ax.plot((cx - scale, cx - scale * 2), (cy + scale, cy + scale), color=color, zorder=3)
            line[0].set_gid(self.uuid + '-' + fragment.name.replace(",", '_') + "-direction")
            patch_dict[fragment.name + "_direction"] = line[0]
            self.data['position'][fragment.name + "_direction"] =\
                (cx - scale, cx - scale * 2), (cy + scale, cy + scale)
            if label:
                ax.text((cx - scale) - 0.32, (cy + scale) + .035, annotation_name)
