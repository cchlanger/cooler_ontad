from setuptools import setup

setup(
    name='use_ontad',
    version='0.1',
    py_modules=['use_ontad'],
    install_requires=[
        'click',
        'bioframe',
        'cooler',
        'numpy',
        'pandas',
    ],
    entry_points='''
        [console_scripts]
        use_ontad=use_ontad:main
    ''',
)
