project('chua', 'c',
  version : '0.1',
  default_options : ['warning_level=3'])

cc = meson.get_compiler('c')

dependencies = [
  # GL usually comes from libglvnd
  dependency('GL'),
  # glu usually comes from mesa
  dependency('glu'),
  # glut usually comes from freeglut
  dependency('glut'),
  # -lm is tricky: https://mesonbuild.com/howtox.html#add-math-library-lm-portably
  cc.find_library('m', required : false),
]

executable(
  'chua',
  sources : [ 'chua.c' ],
  dependencies : dependencies,
  install : true)
