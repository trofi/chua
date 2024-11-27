project('chua', 'c',
  version : '0.1',
  default_options : ['warning_level=3'])

cc = meson.get_compiler('c')

# GL usually comes from libglvnd
gl_dep = dependency('GL')
# glu usually comes from mesa
glu_dep = dependency('glu')
# glut usually comes from freeglut
glut_dep = dependency('glut')

# -lm is tricky: https://mesonbuild.com/howtox.html#add-math-library-lm-portably
m_dep = cc.find_library('m', required : false)

executable(
  'chua',
  sources : [ 'chua.c' ],
  dependencies : [ gl_dep, glu_dep, glut_dep, m_dep ],
  install : true)
