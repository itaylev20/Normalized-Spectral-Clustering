from setuptools import setup, find_packages, Extension

if __name__ == "__main__":
    setup(
        name='myspkmeans', version='0.1.0', author="Itay Lev & Avigail Ben Eliyahu",
        install_requires=['invoke'], packages=find_packages(), license='GPL-2',
        classifiers=[
            'Development Status :: 3 - Alpha',
            'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
            'Natural Language :: English',
            'Programming Language :: Python :: 3 :: Only',
            'Programming Language :: Python :: Implementation :: CPython',
        ],
        ext_modules=[Extension('myspkmeanssp', ['spkmeansmodule.c','spkmeans.c','spk_algorithm.c'])]
    )

