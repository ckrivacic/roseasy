import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

def define_command(module, extras=None):
    entry_point = '{0} = roseasy.commands.{0}:main'.format(module)
    if extras is not None:
        entry_point += ' {0}'.format(extras)
    return entry_point

setuptools.setup(
    name="roseasy-codykrivacic", # Replace with your own username
    version="0.0.1",
    author="Cody Krivacic",
    author_email="krivacic@berkeley.edu",
    description="Rosetta simulations and analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pypa/sampleproject",
    packages=setuptools.find_packages(),
    #packages=find_namespace_packages(include=['roseasy.*']),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    entry_points={
        'console_scripts':[
            'roseasy=roseasy.main:main',
            ],
        'roseasy.commands': [
            define_command('fetch_data'),
            define_command('push_data'),
            define_command('submit'),
            define_command('plot_funnels')
            ],
        }
)
