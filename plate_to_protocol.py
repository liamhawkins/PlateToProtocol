import itertools
from collections import OrderedDict
from opentrons import labware


class PlateFormatError(Exception):
    pass


class PlateLayoutError(Exception):
    pass


class PlateParser:
    def __init__(self, plate_file, format):
        # Make sure plate is set to 96 or 384 well format
        if format not in ['96-flat', 'opentrons-aluminum-block-96-PCR-plate']:
            raise ValueError('Format must be "96-flat" or "opentrons-aluminum-block-96-PCR-plate"')

        self.plate_file = plate_file
        self.format = format
        self.pcr_plate = labware.load(format, '1')
        self.plate_file_list = []
        self.wells = []
        self._rows = []
        self.num_rows = None
        self.num_cols = None
        self._columns = []
        self.pipetting_steps = []

        # Make placeholder wells for specified format
        self._make_plate()

        # Fill plate with data from plate_file
        self._fill_plate()

        # Remove empty rows/columns/wells
        self._clean_plate()

        self.samples = self._get_samples()
        self.targets = self._get_targets()

        self.pipetting_groups = self._get_pipetting_groups()

        self._assign_pipette_group_relations()

        self._assign_pipetting_groups_to_columns()

        self.sample_plate = SamplePlate(self.pipetting_groups)

        self._make_pipetting_steps()



    def __str__(self):
        r = ''
        for row in self.rows():
            r += str(row) + '\n'
        return r

    def rows(self, *args):
        """
        Helper method to access rows of plate by label or index/slice
        """
        row_dict = OrderedDict([(r_label, row) for r_label, row in zip(self.row_labels, self._rows)])
        keys = list(row_dict.keys())

        if not args:
            res = [row_dict[key] for key in keys]
        elif isinstance(args[0], int):
            res = [row_dict[keys[idx]] for idx in args]
        elif isinstance(args[0], str):
            res = [row_dict[idx] for idx in args]
        else:
            raise TypeError

        return res

    def columns(self, *args):
        """
        Helper method to access columns of plate by label or index/slice
        """
        col_dict = OrderedDict([(str(c_label), col) for c_label, col in zip(self.col_labels, self._columns)])
        keys = list(col_dict.keys())

        if not args:
            res = [col_dict[key] for key in keys]
        elif isinstance(args[0], int):
            res = [col_dict[keys[idx]] for idx in args]
        elif isinstance(args[0], str):
            res = [col_dict[idx] for idx in args]
        else:
            raise TypeError
        if len(res) == 1:
            res = res[0]
        return res

    def num_wells(self):
        return len(self.wells)

    def _make_plate(self):
        """
        Create empty placeholder wells, columns, and rows
        """
        if self.format in ['96-flat', 'opentrons-aluminum-block-96-PCR-plate']:
            self.row_labels = list('ABCDEFGH')
            self.col_labels = list(range(1, 12 + 1))
            self.well_vol = 20.0
        else:
            raise NotImplementedError('Format must be "96-flat" or "opentrons-aluminum-block-96-PCR-plate"')

        self.num_rows = len(self.row_labels)
        self.num_cols = len(self.col_labels)

        # Create wells
        for row in self.row_labels:
            for col in self.col_labels:
                self.wells.append(Well(row=row, column=col))
        self.wells = sorted(self.wells, key=lambda x: (x.row, x.column))

        # Populate rows
        for row_label in self.row_labels:
            row_list = [well for well in self.wells if well.row == row_label]
            self._rows.append(Row(row_list))
            # Sort rows on label (A->Z)
            self._rows = sorted(self._rows, key=lambda x: x.label)

        # Populate columns
        for col_label in self.col_labels:
            col_list = [well for well in self.wells if well.column == col_label]
            self._columns.append(Column(col_list))
            # Sort columns on label (1->...)
            self._columns = sorted(self._columns, key=lambda x: x.label)

    def _fill_plate(self):
        """
        Read plate_file and fill in wells with information from file
        """
        # Read in plate file
        with open(self.plate_file, 'r') as f:
            for line in f.read().splitlines():
                self.plate_file_list.append(line.split(','))

        # Find indexes of each column of interest
        self.header = self.plate_file_list[0]
        self.plate_file_wells = self.plate_file_list[1:]
        for i, h in enumerate(self.header):
            if h == 'Row':
                row_index = i
            elif h == 'Column':
                column_index = i
            elif h == 'Sample Type':
                sample_type_index = i
            elif h == 'Replicate #':
                replicate_index = i
            elif h == '*Target Name':
                target_name_index = i
            elif h == '*Sample Name':
                sample_name_index = i
            elif h == 'Starting Quantity':
                starting_quantity_index = i

        # Set proper types to data coming from plate_file then update wells with each entry
        for w in self.plate_file_wells:
            row = w[row_index]
            col = int(w[column_index])
            sample_name = w[sample_name_index] if w[sample_name_index] is not '' else None
            sample_type = w[sample_type_index]
            try:
                replicate_num = int(w[replicate_index])
            except (TypeError, ValueError):
                replicate_num = None
            target_name = w[target_name_index]
            try:
                starting_quantity = float(w[starting_quantity_index])
            except (TypeError, ValueError):
                starting_quantity = None

            # Infer plate format and raise error if plate_file format does not make self.format
            if row not in self.row_labels:
                raise PlateFormatError('Row: {} does not match plate format {}'.format(row, self.format))
            if col not in self.col_labels:
                raise PlateFormatError('Column: {} does not match plate format {}'.format(col, self.format))

            # Update well with corresponding information from plate_file
            self.update_well(row=row,
                             col=col,
                             sample_name=sample_name,
                             sample_type=sample_type,
                             replicate_num=replicate_num,
                             target_name=target_name,
                             starting_quantity=starting_quantity)

        # Check each column to make sure there are not empty wells in the middle of a column.
        # Much harder to support this
        for col in self.columns():
            if col.has_empty_well_in_middle():
                raise PlateLayoutError('Cannot have empty wells in middle of column')

    def _clean_plate(self):
        """
        Remove rows, columns, and wells that are empty
        """
        self._rows = [row for row in self._rows if not row.isempty()]
        self._columns = [col for col in self._columns if not col.isempty()]

        for row in self._rows:
            row.remove_empty_wells()

        for col in self._columns:
            col.remove_empty_wells()

    def update_well(self, row, col, sample_type=None, replicate_num=None, target_name=None, sample_name=None, starting_quantity=None):
        """
        Update information for well corresponding to 'row' and 'col'
        """
        for well in self.wells:
            if well.row == row and well.column == col:
                well.sample_type = sample_type
                well.replicate_num = replicate_num
                well.target_name = target_name
                well.sample_name = sample_name
                well.starting_quantity = starting_quantity

    def _get_samples(self):
        """
        Extract list of samples from wells
        """
        return [Sample(well.sample_name, well.starting_quantity) for well in self.wells]

    def _get_targets(self):
        """
        Extract list of targets from wells
        """
        return [well.target_name for well in self.wells]

    def _get_pipetting_groups(self):
        """
        Creates and returns pipetting groups from non-empty columns. A pipetting group is a series of samples that will
        be pipetted together with the multichannel pipette.
        """
        res = []
        sample_group_dict = OrderedDict()

        # For each non-empty column, create a PipettingGroup containing all non-emtpy wells in that column and add it
        # to or increment its value in sample_group_dict
        for col in self.columns():
            if col.isempty():
                continue
            first_non_empty_well = min(col.non_empty_indexes())
            last_non_empty_well = max(col.non_empty_indexes())
            new_group = PipettingGroup(col[first_non_empty_well:last_non_empty_well+1])
            if new_group not in sample_group_dict:
                sample_group_dict[new_group] = 1
            else:
                sample_group_dict[new_group] += 1

        # Determine which groups are distinct or non-distinct. A distinct group is a PipettingGroup that has no samples
        # in common with any other distinct group. A non-distinct group is a PipettingGroup that is not equal by has at
        # least one sample in common with one of the distinct groups
        distinct_groups, non_distinct_groups = self._get_distinct_and_non_distinct_groups(sample_group_dict)

        # Add distinct groups to return list
        res.extend(distinct_groups)

        # Partition each non-distinct group into new PipettingGroups that are each a subset of one of the
        # distinct_groups. For example if there is a distinct groups containing samples: [1,2,3,4,5,6] and a
        # non-distinct group containing samples: [4,5,6,1,2,3] it would be partitioned into [4,5,6] and [1,2,3]
        for group in non_distinct_groups:
            partitions = self._get_partitions(group, distinct_groups)
            for p in partitions:
                if p not in res:
                    res.append(p)

        return res

    def _get_partitions(self, group, distinct_groups):
        """
        Partition group into PipettingGroups that are subsets of at least one PipettingGroup in distinct_groups
        """
        partition_list = []

        # If group is in distinct_groups return empty list
        if group in distinct_groups:
            return partition_list

        # If group is distinct from all groups in in distinct_groups add to partition_list and return
        if all(sample not in distinct_group for sample in group for distinct_group in distinct_groups):
            partition_list.append(group)
            return partition_list

        # Scan group for longest subset of samples that are in one of the distinct groups
        longest_len = 0
        for first_sample in range(len(group)):
            for last_sample in reversed(range(len(group))):
                candidate_group = group[first_sample:last_sample+1]
                for distinct_group in distinct_groups:
                    if candidate_group in distinct_group and len(candidate_group) > longest_len:
                        longest_len = len(candidate_group)
                        longest_group = candidate_group
                        longest_indexes = (first_sample, last_sample)

        # If partition is not at the beginning of the column, add first sample in column to sample before identified
        # partitions as a new PipettingGroup to partition_list. If this is not a valid PipettingGroup it will be
        # repartitioned below
        if longest_indexes[0] != 0:
            partition_list.append(PipettingGroup(group[0:longest_indexes[0]]))
        # Add identified longest subgroup of samples to partition_list as new PipettingGroup
        partition_list.append(PipettingGroup(longest_group))
        # If identified longest subgroup does not extend to the end of the column, create a new PipettingGroup from
        # remaining samples in the column. If this is not a valid PipettingGroup it will be repartitioned below
        if longest_indexes[1] != len(group)-1:
            partition_list.append(PipettingGroup(group[longest_indexes[1]+1:len(group)]))

        # Scan partition list for invalid PipetteGroups (i.e. a group that is not a subset of one of the
        # distinct_groups) and recursively partition
        repartition_list = []
        for group in partition_list:
            if all([group not in distinct_group for distinct_group in distinct_groups]):
                repartition_list.append(group)
        for group in repartition_list:
            partition_list.remove(group)
            partition_list.extend(self._get_partitions(group, distinct_groups))

        return partition_list

    def _get_distinct_and_non_distinct_groups(self, sample_group_dict):
        """
        Takes dict keyed by sample group, with values of how often that groups occurs in the plate file
        and returns distinct and non-distinct groups. A distinct group is a group that does not have any samples
        in common with the other distinct groups, while a non-distinct group has at least one sample in common with
        one of the distinct groups.
        """
        distinct_groups = []
        non_distinct_groups = []

        # Sort sample_group_dict first by length of the group, then if more than one are longest sort by how
        # many times they occur on the plate
        longest_to_shortest_groups = sorted(sample_group_dict, key=lambda x: (len(x), sample_group_dict[x]))

        # First group is 'distinct' by definition, no other groups have been added so there is nothing
        # for it to be non-distinct from. Subsequent groups are added to distinct_groups if they have no wells
        # in common with any group in distinct_groups, otherwise they are added to non_distinct_groups
        distinct_groups.append(longest_to_shortest_groups[0])
        for group in longest_to_shortest_groups:
            if self._is_distinct_from_groups(group, distinct_groups):
                distinct_groups.append(group)

        for group in longest_to_shortest_groups:
            if group not in distinct_groups:
                non_distinct_groups.append(group)

        return distinct_groups, non_distinct_groups

    def _is_distinct_from_groups(self, group, list_of_groups):
        """
        Check if group is distinct (no samples in common) with all groups in list_of_groups
        """
        for comparator_group in list_of_groups:
            if not self._is_distinct(group, comparator_group):
                return False
        return True

    def _is_distinct(self, group1, group2):
        """
        Check if group1 is distinct (no samples in common) from group2
        """
        samples1 = set(group1.samples)
        samples2 = set(group2.samples)
        if len(samples1.intersection(samples2)) == 0:
            return True
        else:
            return False

    def _assign_pipette_group_relations(self):
        """
        For each PipettingGroup, assign sub- and parent-groups
        """
        for group, comparator_group in [(g1, g2) for g1 in self.pipetting_groups for g2 in self.pipetting_groups if g1 != g2]:
            if comparator_group in group:
                group.add_sub_group(comparator_group)
                comparator_group.add_parent_group(group)

    def _assign_pipetting_groups_to_columns(self):
        """
        Pass pipetting_groups to each col so they can extract their respective PipettingGroups
        """
        for col in self.columns():
            col.match_sample_pipetting_groups(self.pipetting_groups)

    def _make_pipetting_steps(self):
        """
        Iterate over columns and create new PipettingGroup objects with source and destination wells
        """
        for col in self.columns():
            row = ord(col.wells[0].row)
            for pg in col.sample_pipetting_groups:
                new_group = PipettingGroup(pg.samples)
                new_group.source = self.sample_plate.location(pg)
                new_group.dest = self.pcr_plate.wells('{}{}'.format(chr(row), col.label))
                new_group.num_tips = len(pg)
                self.pipetting_steps.append(new_group)
                row += len(pg)


class Well:
    def __init__(self, row, column, sample_type=None, replicate_num=None, target_name=None, sample_name=None, starting_quantity=None):
        self.row = row
        self.column = int(column)
        self.sample_type = sample_type
        try:
            self.replicate_num = int(replicate_num)
        except (ValueError, TypeError):
            self.replicate_num = None
        self.target_name = target_name
        self.sample_name = sample_name
        try:
            self.starting_quantity = float(starting_quantity)
        except (ValueError, TypeError):
            self.starting_quantity = None

    def __str__(self):
        return '({}{})'.format(self.row, self.column)

    def __repr__(self):
        return '{}({}{})'.format(self.__class__.__name__, self.row, self.column)

    def __hash__(self):
        return hash((self.target_name, self.sample_name, self.starting_quantity))

    def __eq__(self, other):
        if isinstance(other, Sample):
            return (self.sample_name == other.sample_name and
                    self.starting_quantity == other.starting_quantity)
        elif isinstance(other, Well):
            return (self.sample_name == other.sample_name and
                    self.target_name == other.target_name and
                    self.starting_quantity == other.starting_quantity)

    def isempty(self):
        return all([x is None for x in [self.sample_type, self.replicate_num, self.target_name, self.sample_name, self.starting_quantity]])


class Sample:
    def __init__(self, sample_name, starting_quantity):
        self.sample_name = sample_name
        self.starting_quantity = starting_quantity

    def __str__(self):
        return '{}:{}-{}'.format(self.__class__.__name__, self.sample_name, self.starting_quantity)

    def __repr__(self):
        return '{}({}, {})'.format(self.__class__.__name__, self.sample_name, self.starting_quantity)

    def __hash__(self):
        return hash((self.sample_name, self.starting_quantity))

    def __eq__(self, other):
        return (self.sample_name == other.sample_name and
                self.starting_quantity == other.starting_quantity)


class PipettingGroup:
    def __init__(self, samples):
        if all(isinstance(r, Sample) for r in samples):
            self.samples = samples
        if all(isinstance(r, Well) for r in samples):
            self.samples = [Sample(r.sample_name, r.starting_quantity) for r in samples]

        self.sub_groups = []
        self.parent_groups = []
        self.num_tips = len(self.samples)
        self.source = None
        self.dest = None
        self.index = 0

    def __eq__(self, other):
        if len(self) != len(other):
            return False
        if isinstance(other, PipettingGroup):
            return all([x == y for x, y in zip(self.samples, other.samples)])
        elif isinstance(other, WellSeries):
            return all([x == y for x, y in zip(self.samples, other.wells)])
        elif all([isinstance(item, (Well, Sample)) for item in other]):
            return all([x == y for x, y in zip(self.samples, other)])
        else:
            raise TypeError("{} is not a PipettingGroup or WellSeries or list of Wells or Samples".format(type(other)))

    def __str__(self):
        s = '{}('.format(self.__class__.__name__)
        for r in self.samples:
            s += str(r) + ', '
        return s[:-2] + ')'

    def __repr__(self):
        sample_str = ''
        for s in self.samples:
            sample_str += repr(s) + ', '
        return '{}([{}])'.format(self.__class__.__name__, sample_str[:-2])

    def __len__(self):
        return len(self.samples)

    def __iter__(self):
        return self

    def __next__(self):
        try:
            result = self.samples[self.index]
        except IndexError:
            self.index = 0
            raise StopIteration
        self.index += 1
        return result

    def __getitem__(self, item):
        return self.samples[item]

    def __hash__(self):
        return hash(tuple(self.samples))

    def __contains__(self, item):
        if isinstance(item, Sample):
            return True if item in self.samples else False

        if isinstance(item, Well):
            return True if Sample(item) in self.samples else False

        # Sliding window comparison
        for i in range(0, len(self)-len(item)+1):
            if self.samples[i:i + len(item)] == item:
                return True
        return False

    def __add__(self, other):
        return PipettingGroup(self.samples + other.samples)

    def add_sub_group(self, group):
        if not isinstance(group, PipettingGroup):
            raise TypeError('Subgroup must be class {}'.format(self.__class__.__name__))
        else:
            self.sub_groups.append(group)

    def add_parent_group(self, group):
        if not isinstance(group, PipettingGroup):
            raise TypeError('Subgroup must be class {}'.format(self.__class__.__name__))
        else:
            self.parent_groups.append(group)

    def has_parents(self):
        return len(self.parent_groups) > 0


class WellSeries:
    def __init__(self, wells=None):
        self.wells = []
        self.sample_pipetting_groups = []
        self.index = 0
        if wells is None:
            self.label = None
        else:
            self.add_wells(wells)

    def __repr__(self):
        well_str = ''
        for w in self.wells:
            well_str += '{}'.format(str(w))
        return "{}({})".format(self.__class__.__name__, well_str)

    def __iter__(self):
        return self

    def __next__(self):
        try:
            result = self.wells[self.index]
        except IndexError:
            self.index = 0
            raise StopIteration
        self.index += 1
        return result

    def __getitem__(self, item):
        return self.wells[item]

    def __len__(self):
        return len(self.wells)

    def __eq__(self, other):
        if len(self) != len(other):
            return False
        if isinstance(other, WellSeries):
            return all([x == y for x, y in zip(self.wells, other.wells)])
        elif isinstance(other, PipettingGroup):
            return all([x == y for x, y in zip(self.wells, other.samples)])
        else:
            raise TypeError("{} is not a PipettingGroup or WellSeries".format(type(other)))

    def __hash__(self):
        return hash(tuple(self.wells))

    def add_wells(self, wells):
        for well in wells:
            self.wells.append(well)

    def isempty(self):
        return all([well.isempty() for well in self])

    def non_empty_indexes(self):
        res = []
        for i, well in enumerate(self):
            if not well.isempty():
                res.append(i)
        return res

    def remove_empty_wells(self):
        self.wells = [well for well in self.wells if not well.isempty()]

    def has_empty_well_in_middle(self):
        empty_index = [index for (index, well) in enumerate(self.wells) if well.isempty()]
        full_index = [index for (index, well) in enumerate(self.wells) if not well.isempty()]

        try:
            min_full = min(full_index)
            max_full = max(full_index)
        except ValueError:
            return False

        if any([index in list(range(min_full, max_full+1)) for index in empty_index]):
            return True
        else:
            return False


class Row(WellSeries):
    def __str__(self):
        s = ''
        for w in self.wells:
            s += str(w)
        return s

    def sort_wells(self):
        self.wells = sorted(self.wells, key=lambda x: x.column)

    def add_wells(self, wells):
        labels = set([x.row for x in wells])
        if len(labels) != 1:
            raise ValueError('Cannot create row from wells with more than 1 row label')
        self.label = labels.pop()

        if self.wells is None:
            self.wells = wells
        else:
            self.wells.extend(wells)
        self.sort_wells()


class Column(WellSeries):
    def __str__(self):
        s = ''
        for w in self.wells:
            s += str(w) + '\n'
        return s

    def sort_wells(self):
        self.wells = sorted(self.wells, key=lambda x: x.row)

    def add_wells(self, wells):
        labels = set([x.column for x in wells])
        if len(labels) != 1:
            raise ValueError('Cannot create row from wells with more than 1 row label')
        self.label = labels.pop()

        if self.wells is None:
            self.wells = wells
        else:
            self.wells.extend(wells)
        self.sort_wells()

    def match_sample_pipetting_groups(self, pipetting_groups):
        for pg in pipetting_groups:
            if pg == self:
                self.sample_pipetting_groups.append(pg)
                return

        # Get all combinations of pipetting_groups whose len is equal to len of column
        combinations = [seq for i in range(len(pipetting_groups), 0, -1) for seq in itertools.permutations(pipetting_groups, i) if sum(len(x) for x in seq) == len(self)]

        # For each combination of PipettingGroups, if combo is equal to column,
        # set that combo as column.pipetting_groups
        for combo in combinations:
            comp_group = combo[0]
            for pg in combo[1:]:
                comp_group += pg
            if comp_group == self:
                self.sample_pipetting_groups.extend(combo)
                return


class SamplePlate:
    def __init__(self, pipetting_groups, plate='96-flat'):
        if plate not in ['96-flat', 'opentrons-aluminum-block-96-PCR-plate']:
            raise NotImplementedError('Only "96-flat" and "opentrons-aluminum-block-96-PCR-plate" are implemented')

        self.pipetting_groups = pipetting_groups
        self.plate = labware.load(plate, '2')
        self.pipetting_group_loc = self._assign_pipetting_group_loc()

    def __str__(self):
        max_len = max([len(str(sample)) for pg in self.pipetting_group_loc.keys() for sample in pg])
        s = ''
        col = 1
        for _ in self.unique_pipetting_groups():
            s += ('  ' + str(col)).ljust(max_len + 1)
            col += 1
        s += '\n'
        row = ord('A')
        for i in range(8):
            s += chr(row + i) + ' '
            for pg in self.unique_pipetting_groups():
                try:
                    s += str(pg[i]).ljust(max_len + 1)
                except IndexError:
                    s += ' '*(max_len + 1)
                    pass
            s += '\n'
        return s

    def _assign_pipetting_group_loc(self):
        pipetting_group_loc = {}
        for pg, col in zip(self.unique_pipetting_groups(), ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12']):
            pipetting_group_loc[pg] = self.plate.wells('A{}'.format(col))

        for pg in self.pipetting_groups:
            if pg.has_parents():
                parent_group = pg.parent_groups[0]
                parent_well_label = str(pipetting_group_loc[parent_group]).replace('<Well ', '').replace('>', '')
                parent_well_row = parent_well_label[0]
                parent_well_col = parent_well_label[1:]
                pg_first_sample = pg[0]

                for index, sample in enumerate(parent_group):
                    if pg_first_sample == sample:
                        pg_row = chr(ord(parent_well_row) + index)
                        pg_col = parent_well_col


                pipetting_group_loc[pg] = self.plate.wells('{}{}'.format(pg_row, pg_col))

        return pipetting_group_loc

    def location(self, pipetting_group):
        return self.pipetting_group_loc[pipetting_group]

    def unique_pipetting_groups(self):
        return [group for group in self.pipetting_groups if not group.has_parents()]



if __name__ == '__main__':
    p = Plate('test_plate4.csv', format='96-flat')

    print("PLATE")
    print(p)

    print('SAMPLE PLATE')
    print(p.sample_plate)

    print('\nPIPETTING STEPS')
    for ps in p.pipetting_steps:
        print('FROM Sample-Plate: {} TO pcr-plate: {} = {}'.format(ps.source, ps.dest, ps))
