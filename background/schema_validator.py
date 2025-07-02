from netCDF4 import Dataset
from enum import Enum, unique, auto
from dataclasses import dataclass
import copy
import numpy as np


@dataclass
class SchemaInfraction:
    @unique
    class InfractionType(Enum):
        DIMENSION_MISSING = auto()
        DIMENSION_WRONG_SIZE = auto()
        VARIABLE_MISSING = auto()
        VARIABLE_WRONG_DTYPE = auto()
        VARIABLE_WRONG_DIMENSIONS = auto()
        ATTRIBUTE_MISSING = auto()
        ATTRIBUTE_WRONG_VALUE = auto()
        GROUP_MISSING = auto()

    infraction_type: InfractionType  # What kind of error was found?
    infraction_element: str  # Name of group, variable, or dimension in question.
    infraction_description: str  # Long form description of the error.


def _deep_equals(a, b):
    """
    When working with NetCDF attributes, they can either be single objects, or they can be lists
    or dictionaries. We want to be able to compare regardless, so pull out this comparison here.

    Input:
        a: first object
        b: second object

    Output: Boolean, true if the two objects are equal, false otherwise.

    It really seems like this should be built into Python or a library or something - if it is,
    I can't find it. Specifically, we need support for Numpy conversion like this - NetCDF will
    read the list as a Numpy array, while Yaml reads as a plain old list.
    """
    # Caveat, if a or b is a Numpy array, force to a list.
    if isinstance(a, np.ndarray):
        a = a.tolist()
    if isinstance(b, np.ndarray):
        b = b.tolist()
    # A is a List.
    if isinstance(a, list):
        # If B is not a list, shortcut, they can't be equal.
        # Caveat: Numeric li
        if not isinstance(b, list):
            print("B not list")
            return False
        # Check lists are the same length
        if len(a) != len(b):
            return False
        # If any elements are not equal (order matters), return False.
        for a_elem, b_elem in zip(a, b):
            if not _deep_equals(a_elem, b_elem):
                return False
        # Two lists with equal elements are equal.
        return True
    # A is a dictionary.
    elif isinstance(a, dict):
        # If B is not a dict, shortcut, they can't be equal.
        if not isinstance(b, dict):
            return False
        # Check dictionaries are the same length
        if len(a) != len(b):
            return False
        # Check that all of A's elements are in B.
        for key in a.keys():
            if not (key in b) or not _deep_equals(a[key], b[key]):
                return False
        return True
    # A is a floating point number.
    elif isinstance(a, (float, np.floating)):
        # If B is not a float type, shortcut, they can't be equal.
        if not isinstance(b, (float, np.floating)):
            return False
        # Check for equality within a tolerance.
        return np.allclose(a, b)
    # If not a dict, list, or float, compare directly.
    return a == b


def validate_against_schema(dataset: Dataset, schema: dict) -> tuple[bool, list[SchemaInfraction]]:
    """
    Validate the given NetCDF4 dataset against a schema. See msat_netcdf/README.md for full details
    on schema structure.  A dataset is considered breaking the schema if it does not contain the
    necessary data, but additional data is allowed (but should not be relied on)

    Input:
        dataset: The NetCDF4 dataset to validate.
        schema: A loaded YAML describing the schema.

    Output:
        True, []: if the dataset fits the schema
        False, list_of_infractions: if the dataset does not fit the schema.
    """
    infractions = []
    # Avoid modifying the schema as we pop elements from it.
    schema = copy.deepcopy(schema)
    # Verify dimensions in this group. Actual dataset may have more dimensions.
    expected_dimensions = schema.pop("dimensions", {})
    for dim, size in expected_dimensions.items():
        if dim not in dataset.dimensions:
            infractions.append(
                SchemaInfraction(
                    SchemaInfraction.InfractionType.DIMENSION_MISSING,
                    dim,
                    f"Dataset expected to have dimension {dim} in group {dataset.name}.",
                )
            )
        # size = 0 is a signal that any size is acceptable.
        # This follows netCDF4 convention for unlimited dimensions.
        elif size != 0 and size != dataset.dimensions[dim].size:
            infractions.append(
                SchemaInfraction(
                    SchemaInfraction.InfractionType.DIMENSION_WRONG_SIZE,
                    dim,
                    f"Dimension {dim} expected to have size {size} but has size"
                    f" {dataset.dimensions[dim].size}.",
                )
            )

    # Read the schema for variables, groups, or attributes
    variables = []
    groups = []
    attributes = []
    for key, value in schema.items():
        if isinstance(value, dict):  # variable or group
            if "type" in value:
                variables.append((key, value))
            else:
                groups.append((key, value))
        else:
            attributes.append((key, value))

    # Verify existence of all variables.
    for name, metadata in variables:
        # Check if variable is present at all.
        if name not in dataset.variables:
            infractions.append(
                SchemaInfraction(
                    SchemaInfraction.InfractionType.VARIABLE_MISSING,
                    name,
                    f"Dataset expected to have variable {name} in group {dataset.name}.",
                )
            )
            continue  # Variable missing, no point checking other attributes.
        actual_variable = dataset[name]

        # Pop special fields before verifying others as attributes.
        # Verify datatype
        expected_datatype = metadata.pop("type")
        if not np.issubdtype(actual_variable.dtype, expected_datatype):
            infractions.append(
                SchemaInfraction(
                    SchemaInfraction.InfractionType.VARIABLE_WRONG_DTYPE,
                    name,
                    f"Variable {name} expected to have dtype {expected_datatype} but has type"
                    f" {actual_variable.dtype}.",
                )
            )
        # Verify dimensions and order
        expected_dimensions = tuple(metadata.pop("dimensions", ()))
        if actual_variable.dimensions != expected_dimensions:
            infractions.append(
                SchemaInfraction(
                    SchemaInfraction.InfractionType.VARIABLE_WRONG_DIMENSIONS,
                    name,
                    f"Variable {name} expected to have dimensions {expected_dimensions} but has"
                    f" dimensions {actual_variable.dimensions}.",
                )
            )
        # Ignore DictionaryKey as it's not expected in the dataset.
        metadata.pop("dictionary_key", None)
        for attribute, value in metadata.items():
            if attribute not in actual_variable.ncattrs():
                infractions.append(
                    SchemaInfraction(
                        SchemaInfraction.InfractionType.ATTRIBUTE_MISSING,
                        name,
                        f"Variable {name} expected to have attribute {attribute}.",
                    )
                )
            # We don't want to trigger on blank values for attributes.
            elif value and not _deep_equals(value, actual_variable.getncattr(attribute)):
                infractions.append(
                    SchemaInfraction(
                        SchemaInfraction.InfractionType.ATTRIBUTE_WRONG_VALUE,
                        name,
                        f"Attribute {attribute} on {name} expected to have value {value}, but is"
                        f" {actual_variable.getncattr(attribute)}.",
                    )
                )
    # Now check Attributes on the group
    for attribute, value in attributes:
        if attribute not in dataset.ncattrs():
            infractions.append(
                SchemaInfraction(
                    SchemaInfraction.InfractionType.ATTRIBUTE_MISSING,
                    name,
                    f"Variable {name} expected to have attribute {attribute}.",
                )
            )
        elif value and not _deep_equals(value, dataset.getncattr(attribute)):
            infractions.append(
                SchemaInfraction(
                    SchemaInfraction.InfractionType.ATTRIBUTE_WRONG_VALUE,
                    name,
                    f"Attribute {attribute} on {name} expected to have value {value}, but is"
                    f" {dataset.getncattr(attribute)}.",
                )
            )

    # Now, check groups and recurse into them, adding all their infractions too.
    for name, value in groups:
        if name not in dataset.groups:
            infractions.append(
                SchemaInfraction(
                    SchemaInfraction.InfractionType.GROUP_MISSING,
                    name,
                    f"Group {dataset.name} expected to have subgroup {name}.",
                )
            )
        else:
            _, subinfractions = validate_against_schema(dataset[name], value)
            infractions += subinfractions

    return len(infractions) == 0, infractions
