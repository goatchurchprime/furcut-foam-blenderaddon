bl_info = {
        'name': 'Furcut foadm',
        'author': 'Julian Todd',
        'version': (1, 0),
        'blender': (5, 0, 0),
        'category': 'CNC',
        'description': 'CNC foam cutting',
        'location': 'Object > Seams to Sewing Pattern > ...',
        'wiki_url': 'https://blenderartists.org/t/1248713'}

from . import op_upcut
from . import op_upslice
from . import op_post
if "bpy" in locals():
    import importlib
    importlib.reload(op_test)
    importlib.reload(op_upcut)
    importlib.reload(op_upslice)
    importlib.reload(op_post)
else:
    from . import op_test
    from . import op_upcut
    from . import op_upslice
    from . import op_post

# need to edit the __init__.py file to trigger reload when addin by selecting "seams to sewing patterns" from edit->preferences and then disable and re-enable
print("RELOADING FURCUT ADDON *** YdfaddFssfds-ddsddgdddddsssssd-sddsdeds-sshdsgsdsdsdsssdassfY ***")

import bpy
from bpy.types import Menu

def menu_func(self, context):
    lay_out = self.layout
    lay_out.operator_context = 'INVOKE_REGION_WIN'
    lay_out.separator()
    lay_out.menu("VIEW3D_MT_object_furcut_menu", text="Furcut Foam")
    
class VIEW3D_MT_object_furcut_menu(Menu):
    bl_idname = "VIEW3D_MT_object_furcut_menu"
    bl_label = "Furcut FFoam"

    def draw(self, context):
        layout = self.layout
        layout.operator("object.furcut_test", text="Furcut foam test", icon="OUTLINER_DATA_SURFACE")
        layout.operator("object.furcut_upcut", text="Furcut foam upcut", icon="OUTLINER_DATA_SURFACE")
        layout.operator("object.furcut_upslice", text="Furcut foam upslice", icon="OUTLINER_DATA_SURFACE")
        layout.operator("object.furcut_post", text="Furcut foam post", icon="OUTLINER_DATA_SURFACE")
        layout.separator()

classes = [
    VIEW3D_MT_object_furcut_menu,
    op_test.TestTest,
    op_upslice.Upslice,
    op_upcut.Upcut,
    op_post.PostProc,
    ]

def register():
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)
    bpy.types.VIEW3D_MT_object.append(menu_func)

def unregister():
    bpy.types.VIEW3D_MT_object.remove(menu_func)
    from bpy.utils import unregister_class
    for cls in reversed(classes):
        unregister_class(cls)

if __name__ == "__main__":
    register()
