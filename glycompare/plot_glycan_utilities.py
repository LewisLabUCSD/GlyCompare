import matplotlib.pyplot as plt

from glypy.plot.draw_tree import *
from glypy.plot.buchheim import buchheim
from glypy.plot.topological_layout import layout as topological
from glypy.plot import cfg_symbols, iupac_symbols



# def

def _bounding_box_to_limits(dtree):
    scale_up = mtransforms.Affine2D().scale(1.1).get_matrix()
    points = dtree._compute_bounding_box()
    stretched = []
    for x, y in points:
        x, y, w = scale_up.dot(np.array((x, y, 1)))
        x /= w
        y /= w
        stretched.append((x, y))
    points = stretched
    xs, ys = zip(*points)
    min_x, max_x = min(xs) - 2, max(xs) + 2
    min_y, max_y = min(ys) - 2, max(ys) + 2
    return np.array((min_x, max_x)), np.array((min_y, max_y))


line_to = cfg_symbols.CFGNomenclature().line_to

nomenclature_map = {
    "cfg": cfg_symbols.CFGNomenclature(),
    'iupac': iupac_symbols.IUPACTextSymbolicNomenclature(),
    "snfg": snfg_symbols.SNFGNomenclature()
}

layout_map = {
    "balanced": BalancedTreeLayout,
    "topological": TopologicalTreeLayout
}

# DEFAULT_SYMBOL_SCALE_FACTOR = 0.16

#: :data:`special_cases` contains a list of names for
#: special case monosaccharides
_special_case_fucose = monosaccharides.Fuc
_special_case_fucose.configuration = None
_special_case_fucose.anomer = None
_special_case_xylose = monosaccharides.Xyl
_special_case_xylose.configuration = None
_special_case_xylose.anomer = None

special_cases = [_special_case_fucose, _special_case_xylose]

anomer_symbol_map = {
    constants.Anomer.alpha: r'\alpha',
    constants.Anomer.beta: r'\beta',
    constants.Anomer.uncyclized: r'o',
    constants.Anomer.x: r'?'
}



def plot_glycan_profile(a_profile, glycan_dict):
    _a = len(a_profile.keys())
    _len = 5
    _r = divmod(_a, 5)[0] + 1 if divmod(_a, 5)[1] != 0 else (divmod(_a, 5)[0])
    fig, axes = plt.subplots(_r, 5, squeeze=False)

    fig_size = (divmod(_a, 5)[0] + 1) * 3 if divmod(_a, 5)[1] != 0 else (divmod(_a, 5)[0]) * 3
    fig.set_size_inches(16, fig_size)
    #     print(axes)
    #     fig.set_size_inches(16,5)
    _count = 0
    for i in sorted(a_profile.keys()):
        #             print(i)
        #             print(divmod(_count, _a))
        _x, _y = divmod(_count, _len)
        plot_glycan(glycan_dict[a_profile[i]], ax=axes[_x][_y], center=True, title=str(i))
        _count += 1


def plot_glycan_list(glycan_obj_list, idex_list=[], title='Glycans', addr=''):
    """
    :param glycan_obj_list: a list of Glycan
    :param idex_list: a list of index for glycan (i.e. motif indexes)
    :param title: str title
    :return:
    """
    _a = len(glycan_obj_list)
    _len = 5
    _r = divmod(_a, 5)[0] + 1 if divmod(_a, 5)[1] != 0 else (divmod(_a, 5)[0])
    fig, axes = plt.subplots(_r, 5, squeeze=False)
    plt.title(title)
    fig_size = (divmod(_a, 5)[0] + 1) * 2.5 if divmod(_a, 5)[1] != 0 else (divmod(_a, 5)[0]) * 2.5
    fig.set_size_inches(15, fig_size)
    #     print(axes)
    #     fig.set_size_inches(16,5)
    _count = 0
    for idex, i in enumerate(glycan_obj_list):
        #             print(i)
        #             print(divmod(_count, _a))
        _x, _y = divmod(_count, _len)
        if not idex_list:
            plot_glycan(i, ax=axes[_x][_y], center=True, title=str(idex))
        else:
            plot_glycan(i, ax=axes[_x][_y], center=True, title=str(idex_list[idex]))
        _count += 1
    if addr == '':
        pass
    else:
        plt.savefig(addr)


def plot_glycan_dict(a_glycan_dict, single=True):

    a_keys = list(a_glycan_dict.keys())
    if single:
        for i in a_keys:
            plot_glycan(a_glycan_dict[i],title=i)
    else:
        plot_glycan_list([a_glycan_dict[i] for i in a_keys], idex_list=a_keys)


def plot_glycan(tree, title='', addr='', at=(0, 0), ax=None, orientation='h', center=True, label=False,
                symbol_nomenclature='cfg', layout='balanced', layout_args=None, **kwargs):
    '''
    Draw the parent outlink position and the child anomer symbol
    Parameters
    ----------
    tree: Glycan or Monosaccharide
        The saccharide structure to draw.
    orientation: str
        A string corresponding to `h` or `horizontal` will draw the glycan horizontally
        from right to left. `v` or `vertical` will draw the glycan from bottom to top.
        Defaults to `h`
    at: tuple
        The x, y coordinates at which to draw the glycan's reducing_end. Defaults to `(0, 0)`
    ax: :class:`matplotlib.axes.Axis`
        A matplotlib axis on which to draw. Useful for adding glycans to existing figures. If
        an axis is not provided, a new figure will be created and used.
    center: bool
        Should the plot limits be centered around this glycan? Defaults to |False| but will be
        be set to |True| if `ax` is |None|
    label: bool
        Should the bond annotations for `tree` be drawn? Defaults to |False|
    scale: (float, float) or float
        Scale coefficients. Pass a pair to scale x and y dimensions respectively.
    '''

    scale = DEFAULT_SYMBOL_SCALE_FACTOR
    transform_scale = kwargs.pop("scale", (0.8,))
    try:
        transform_scale = tuple(transform_scale)
    except Exception:
        transform_scale = (transform_scale, transform_scale)

    if isinstance(symbol_nomenclature, basestring):
        symbol_nomenclature = nomenclature_map.get(symbol_nomenclature)
    elif symbol_nomenclature is None:
        symbol_nomenclature = snfg_symbols.SNFGNomenclature()

    if isinstance(layout, basestring):
        layout = layout_map.get(layout)

    tree_root = root(tree)
    dtree = DrawTreeNode(tree_root)

    layout_algorithm = layout(dtree)
    layout_algorithm.layout(**(layout_args or {}))

    fig = None
    # Create a figure if no axes are provided
    if ax is None:
        fig, ax = plt.subplots()
        # at = (0, 0)
        center = True
        ax.axis('off')
    (layout_transform,
     substituent_transform,
     _) = symbol_nomenclature.get_layout_transform(orientation)
    layout_algorithm.transform(layout_transform)
    dtree.draw(at=at, ax=ax, scale=scale, label=False,
               symbol_nomenclature=symbol_nomenclature,
               annotation_transform=substituent_transform,
               **kwargs)
    dtree.axes = ax
    dtree.data['orientation'] = orientation
    dtree.set_transform(mtransforms.Affine2D())
    if label:
        dtree.draw_linkage_annotations(
            at=at, ax=ax, scale=scale,
            symbol_nomenclature=symbol_nomenclature,
            **kwargs)
    if orientation in {"h", "horizontal"}:
        pass
        # dtree.set_transform(mtransforms.Affine2D().rotate_deg(90).scale(*transform_scale))
        # dtree.update_text_position(-90)
    # If the figure is stand-alone, center it
    if fig is not None or center:
        if title != "":
            ax.title.set_text(title)
        xlim, ylim = _bounding_box_to_limits(dtree)
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_facecolor('white')

        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
    if addr == '':
        pass
    else:
        plt.savefig(addr, transparent=True, format='eps', dpi='figure')
    return (dtree, ax)


#
# if __name__ == '__main__':
#     output_glycan_motif_vec_to_file()


if __name__ == '__main__':
    output_glycan_motif_vec_to_file()
