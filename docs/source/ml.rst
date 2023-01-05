``lipidoz.ml``
=======================================
This subpackage contains utilities for performing double bond identification using machine learning.


Module Reference
-------------------------------------------------------------

``lipidoz.ml.data``
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.. autofunction:: lipidoz.ml.data.load_preml_data

.. autofunction:: lipidoz.ml.data.load_ml_targets

.. autofunction:: lipidoz.ml.data.split_true_and_false_preml_data

.. autofunction:: lipidoz.ml.data.preml_to_ml_data

.. autofunction:: lipidoz.ml.data.get_dataloaders_for_ml


``lipidoz.ml.models``
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ML Model Base Class Methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: lipidoz.ml.models._Model.load

.. autofunction:: lipidoz.ml.models._Model.save

.. autofunction:: lipidoz.ml.models._Model.train

.. autofunction:: lipidoz.ml.models._Model.predict

.. autofunction:: lipidoz.ml.models._Model.predict_proba

ResNet18 model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: lipidoz.ml.models.resnet18.ResNet18

.. autofunction:: lipidoz.ml.models.resnet18.ResNet18.__init__


``lipidoz.ml.view``
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.. autofunction:: lipidoz.ml.view.plot_preml_example

.. autofunction:: lipidoz.ml.view.plot_ml_example
