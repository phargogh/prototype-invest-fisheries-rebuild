import logging

import numpy
import pandas
from natcap.invest import datastack

logging.basicConfig(level=logging.INFO)
LOGGER = logging.getLogger(__name__)


def main():
    # workspace dir
    # populations params file
    #    * same structure as the fisheries population parameters table.
    # whether population classes are sex-specific
    #    * boolean
    # Habitat dependency parameters file
    #    * Columns: age/stage names
    #    * Rows: habitat names
    #    * Cell values: habitat dependencies (float, 0-1)
    # habitat area change file
    #    * Columns: subregion names
    #    * Rows: habitat names
    #    * Cell values: % changes in habitat area by subregion.
    # gamma (numeric parameter)
    #    * Float value between 0 and 1
    #    * Describes the relationship between the change in habitat area and a
    #      change in survival of age/stage dependency on that habitat.



if __name__ == '__main__':
    main()


