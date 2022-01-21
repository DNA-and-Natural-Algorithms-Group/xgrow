cc = meson.get_compiler('c')
m_dep = cc.find_library('m', required : false)

x11 = dependency('x11')

pymod = import('python')
python = pymod.find_installation('python3')
python_dir = python.get_install_dir()

xgrow = executable('_xgrow', '../src/xgrow.c', '../src/grow.c',
           dependencies: [x11, m_dep], install: true, install_dir: python_dir + '/xgrow')

library('xgrow', '../src/grow.c', dependencies: m_dep)



python_sources = [
  'xgrow.py',
  'stxg.py',
  'tileset.py',
  'xcolors.py',
]

python.install_sources(python_sources,
  pure : false,
  subdir : 'xgrow'
)
