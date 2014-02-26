.. _Demos:

Demos
==================================
To get started, we recommend starting with the demos. To get access to all the
demos, execute the following command in a terminal window:

.. code-block:: bash

   cbcflow-get-demos
   
To list and run all the demos, execute

.. code-block:: bash

   cd cbcflow-demos/demo
   ./cbcflow_demos.py --list
   ./cbcflow_demos.py --run
   
   
If you have downloaded the development version, it is sufficient to download the
demo data in the root folder of the repository:

.. code-block:: bash

   cbcflow-get-data

If you are
unfamiliar with FEniCS, please refer to the `FEniCS Tutorial <http://fenicsproject.org/documentation/tutorial/>`_
for the FEniCS-specifics of these demos.


**Documented demos**:

.. toctree::
    :maxdepth: 1
    
    FlowAroundCylinder/documentation
    Womersley3D/documentation
    Replay/documentation
    Restart/documentation
    
