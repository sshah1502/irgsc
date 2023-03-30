import pdoc
modules = ['src.irgsctool.__init__', 'src.irgsctool._get_data',\
    'src.irgsctool._read_data', 'src.irgsctool._sam', 'src.irgsctool._sgc',\
    'src.irgsctool._extinction_correction', 'src.irgsctool._fitting', 'src.irgsctool._validate']  # Public submodules are auto-imported
context = pdoc.Context()
modules = [pdoc.Module(mod, context=context)
           for mod in modules]
pdoc.link_inheritance(context)
def recursive_htmls(mod):
    yield mod.name, mod.html()
    for submod in mod.submodules():
        yield from recursive_htmls(submod)
for mod in modules:
    for module_name, html in recursive_htmls(mod):
        ...  # Process