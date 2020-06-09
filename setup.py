from setuptools import setup

setup(
    name='cooler_ontad',
    version='0.1',
    py_modules=['cooler_ontad'],
    install_requires=[
        'click',
        'bioframe',
        'cooler',
        'numpy',
        'pandas',
    ],
    entry_points='''
        [console_scripts]
        cooler_ontad=cooler_ontad:main
    ''',
)
