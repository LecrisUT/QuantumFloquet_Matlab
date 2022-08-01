#!/usr/bin/env python
import os
import re
import warnings
import argparse

import yaml
import xml.etree.ElementTree as ET
from jinja2 import Environment, FileSystemLoader

currYml = {}
fullYml = {}
currKey = ''
currCls = ''
currPkg = []
currPath = ''
currDocPath = ''
rootPath = ''
currIndex = {}
env = Environment(loader=FileSystemLoader('.'),
                  trim_blocks=True, lstrip_blocks=True)
docroot = "../ref"
htmlroot = "../html"


# Helpers
def existYAML(value, onlyRelative=False, *relPath):
    global currPath, rootPath, currPkg, currCls
    value = value.split('.')
    pkg = currPkg
    cls = currCls
    if len(value) == 1:
        # Checking YAML with respect to current requester
        value = value[0]
    else:
        # Checking YAML of a FQN
        # Check if constructor or package.function
        if value[-1] == value[-2]:
            value.pop()
            if existYAML('.'.join(value), *relPath):
                # Is constructor
                cls = value[-1]
                pkg = value[:-1]
                value = 'constructor'
            else:
                # Is package.function
                pkg = value
                cls = ''
                value = value[-1]
        else:
            pkg = value[:-1]
            cls = ''
            value = value[-1]
    # TODO: Special case of cls, check if inside cls.yaml file
    if onlyRelative:
        paths = [os.path.join(currPath, *relPath, value + '.yaml'),
                 os.path.join(currPath, *relPath, value + '.yml')]
    else:
        # Prepend appropriate prefixes to path search
        pkg = ['+' + p for p in pkg]
        if cls:
            cls = '@' + cls
        paths = [os.path.join(currPath, *relPath, value + '.yaml'),
                 os.path.join(currPath, *relPath, value + '.yml'),
                 os.path.join(rootPath, *pkg, cls, value + '.yaml'),
                 os.path.join(rootPath, *pkg, cls, value + '.yml')]
    if any(os.path.exists(testPath := a) for a in paths):
        return True, testPath
    return False, ""


def existPackage(value, onlyRelative=False, onlyAbsolute=False):
    global currPath, rootPath
    value = value.split('.')
    # Prepend appropriate prefixes to path search
    value = ['+' + v for v in value]
    if onlyAbsolute:
        paths = [os.path.join(rootPath, *value)]
    elif onlyRelative:
        paths = [os.path.join(currPath, *value)]
    else:
        paths = [os.path.join(currPath, *value),
                 os.path.join(rootPath, *value)]
    if any(os.path.exists(testPath := a) for a in paths):
        return True, testPath
    return False, ""


# Custom filters
def isHTML(value):
    # Regex: </?\s*[a-z-][^>]*\s*>|(\&(?:[\w\d]+|#\d+|#x[a-f\d]+);)
    # https://regex101.com/r/Va987X/1
    # https://stackoverflow.com/a/70585604
    return re.match(r"</?\s*[a-z-][^>]*\s*>|(\&(?:[\w\d]+|#\d+|#x[a-f\d]+);)", value)


def asOutputs(values, withLinks=False):
    if not values:
        return ""
    if type(values) is str:
        # Special case:
        if values == 'all':
            return "[___] = "
        values = [values]
    if type(values) is list:
        if not all(type(val) is str for val in values):
            raise TypeError("Input type is not string")
        if len(values) == 1:
            if withLinks:
                return f"<a class=\"intrnllnk\" href=\"#out-{values[0]}\">{values[0]}</a> = "
            return f"{values[0]} = "
        out = "["
        for val in values:
            if withLinks:
                out += f"<a class=\"intrnllnk\" href=\"#out-{val}\">{val}</a>,"
            else:
                out += f"{val},"
        out = out.rstrip(',')
        out += '] = '
        return out
    else:
        raise TypeError("Input type unrecognized")


def asInputs(values, withLinks=False):
    global currYml
    if not values:
        return ""
    if type(values) is str:
        # Special case:
        if values == 'all':
            return "(___)"
        values = [values]
    if type(values) is list:
        if not all(type(val) is str for val in values):
            raise TypeError("Input type is not string")
        out = "("
        repeating = []
        nameValued = False
        for val in values:
            # Special case for nameValued
            val = val.split('=')
            if len(val) == 2:
                nameValued = True
                if not currYml['nameValue'].get(val[0]):
                    raise ValueError(f"Input {val[0]} is not a recognized name-value input")
                if withLinks:
                    out += f"<a class=\"intrnllnk\" href=\"#nameValue-{val[0]}\">{val[0]}</a>={val[1]},"
                else:
                    out += f"{val[0]}={val[1]},"
            elif len(val) == 1:
                val = val[0]
                # Check if it is not a regular inputs
                if not currYml['inputs'].get(val):
                    if not currYml['nameValue'].get(val):
                        raise ValueError(f"Input {val} is neither a regular input nor a name-value input")
                    if withLinks:
                        out += f"<a class=\"intrnllnk\" href=\"#nameValue-{val}\">{val}</a>={val},"
                    else:
                        out += f"{val}={val},"
                else:
                    if nameValued:
                        raise SyntaxError("Regular or repeating inputs cannot be after name-value inputs")
                    if currYml['inputs'][val].get('repeating'):
                        repeating += [val]
                    else:
                        if repeating:
                            raise SyntaxError("Inputs are ill-structured. Repeating inputs have to be at the end, before name-value inputs")
                        if withLinks:
                            out += f"<a class=\"intrnllnk\" href=\"#inp-{val}\">{val}</a>,"
                        else:
                            out += f"{val},"
            else:
                raise ValueError(f"Input {'='.join(val)} is ill-defined")
        if repeating:
            for val in repeating:
                if withLinks:
                    out += f"<a class=\"intrnllnk\" href=\"#inp-{val}\">{val}1</a>,"
                else:
                    out += f"{val}1,"
            out += "...,"
            for val in repeating:
                if withLinks:
                    out += f"<a class=\"intrnllnk\" href=\"#inp-{val}\">{val}N</a>,"
                else:
                    out += f"{val}N,"
        out = out.rstrip(',')
        out += ')'
        return out
    else:
        raise TypeError("Input type unrecognized")


def getDataType(value):
    if value.get('dataType'):
        return value['dataType']
    out = ""
    if value.get('class'):
        out += f"{value['class']} "
    if value.get('dims'):
        dims = value['dims']
        if type(dims) is not list:
            dims = [dims]
        for ind, ndims in enumerate(dims):
            if type(ndims) is not int or ndims < 0:
                raise NameError("Input dims is not a non-negative integer")
            if ndims == 0:
                out += "scalar "
            elif ndims == 1:
                out += "vector "
            elif ndims == 2:
                out += "matrix "
            else:
                out += f"{ndims}D-array "
            if ind != len(dims) - 1:
                out += "or "
    if not out:
        out = "arbitrary "
    return out


def getTypes(value):
    types = []
    if value.get('class'):
        types += [value['class']]
    if value.get('dims'):
        dims = value['dims']
        if type(dims) is list:
            dims = max(dims)
        if type(dims) is not int or dims < 0:
            raise NameError("Input dims is not a non-negative integer")
        if dims == 0:
            types += ["scalar"]
        elif dims == 1:
            types += ["vector"]
        elif dims == 2:
            types += ["2d"]
        elif dims == 3:
            types += ["3d"]
    types = list(dict.fromkeys(types))
    return types


def getFQN(value):
    global currPkg, currKey
    if value:
        if value == 'constructor':
            # TODO: Currently assuming currKey is correct constructor format
            return '.'.join([*currPkg, currKey])
        else:
            return '.'.join([*currPkg, value])
    return '.'.join([*currPkg, currKey])


def getPurpose(value):
    global currYml, fullYml
    # Check if purpose is in current yaml
    if type(currYml['methods'][value]) is dict and \
            currYml['methods'][value].get('purpose') and \
            type(currYml['methods'][value]['purpose']) is str:
        return currYml['methods'][value]['purpose']
    # Search other yamls
    if (testPath := a for _, a in existYAML(value)):
        with open(testPath) as fl:
            yml = yaml.load(fl, yaml.Loader)
        yml = list(yml.values())[0]
        if yml.get('purpose') and type(yml['purpose']) is str:
            return yml['purpose']
        else:
            warnings.warn(f"No purpose field was found in {testPath}")
    # Special case for methods in class template
    elif currCls:
        value = value.split('.')
        if len(value) == 1:
            value = value[0]
            if value == 'constructor' and fullYml.get(f"{currCls}.{currCls}"):
                yml = fullYml[f"{currCls}.{currCls}"]
                if yml.get('purpose') and type(yml['purpose']) is str:
                    return yml['purpose']
                else:
                    warnings.warn(f"No purpose field was found in {currCls}.{currCls}")
            elif fullYml.get(f"{currCls}.{value}"):
                yml = fullYml[f"{currCls}.{value}"]
                if yml.get('purpose') and type(yml['purpose']) is str:
                    return yml['purpose']
                else:
                    warnings.warn(f"No purpose field was found in {currCls}.{value}")
        elif len(value) == 2:
            if currCls != value[0]:
                raise ValueError(f"Referencing different classes ({currCls} =/= {value[0]}) for method {value[1]}")
            value = value[1]
            if fullYml.get(f"{currCls}.{value}"):
                yml = fullYml[f"{currCls}.{value}"]
                if yml.get('purpose') and type(yml['purpose']) is str:
                    return yml['purpose']
                else:
                    warnings.warn(f"No purpose field was found in {currCls}.{value}")
        else:
            warnings.warn(f"Trying to reference external object, but yaml is not found")
    return ""


def getPackage():
    global currPkg
    return '.'.join(currPkg)


def getHRef(value):
    global currPkg, currKey, currYml, currPath, rootPath, currDocPath
    # Find the current documentation path
    depth = len(currDocPath.split(os.sep))
    # These documentation references are under ref/... so count it as well
    # TODO: Do not hardcode
    depth = depth + 1
    # Check for special cases
    if value.startswith('http'):
        # External links return as is
        return value
    if value.endswith('.html'):
        # Local html link are relative current documentation root
        return '/'.join(['..'] * depth) + '/' + value
    # Otherwise: main logic
    # Check if it is a local documentation
    value = value.split('/')
    if len(value) == 1:
        subref = ""
    elif len(value) == 2:
        subref = value.pop()
    else:
        raise ValueError(f"Reference: {'/'.join(value)} is ill-formatted")
    value = value[0]
    value = value.split('.')
    if len(value) == 1:
        value = value[0]
        if existYAML(value, onlyRelative=True):
            # If class, function, method, etc. get the current link
            ref = f"{value}.html"
        elif existPackage(value, onlyRelative=True):
            # If package
            ref = f"{value}/index.html"
        else:
            # Otherwise pass to matlab:doc
            return f"matlab:doc(&#34;{value}&#34;)"
    else:
        if existYAML('.'.join(value)):
            ref = '/'.join(['..'] * (depth - 1)) + '/' + '/'.join(value) + '.html'
        elif existPackage('.'.join(value), onlyAbsolute=True):
            ref = '/'.join(['..'] * (depth - 1)) + '/' + '/'.join(value) + '/index.html'
        else:
            return f"matlab:doc(&#34;{'.'.join(value)}&#34;)"
    # TODO: Resolve reference if input, output, etc.
    # Default pass reference as is
    if subref:
        subref = '#' + subref
    return ref + subref


# Handle Template
def checkFunction():
    global currYml, currPkg, currKey
    valid = True
    key = '.'.join([*currPkg, currKey])
    if type(currYml.get('syntax')) is not dict:
        warnings.warn(f"{key} lacks field syntax or is not object")
        valid = False
    elif type(currYml['syntax'].get('groups')) is not list or \
            not all(type(d) is dict for d in currYml['syntax']['groups']):
        warnings.warn(f"{key} lacks field syntax.groups or is not list of objects")
        valid = False
    elif not all(type(d.get('syntax')) is list and
                 all(type(s) is dict for s in d['syntax'])
                 for d in currYml['syntax']['groups']):
        warnings.warn(f"{key} syntax.groups[].syntax fields is missing or not an object")
        valid = False
    if currYml.get('inputs'):
        if type(currYml['inputs']) is not dict:
            warnings.warn(f"{key} inputs is not an object")
            valid = False
        elif not all(d.get('purpose') and type(d['purpose']) is str
                     for d in currYml['inputs'].values()):
            warnings.warn(f"{key} input object is missing purpose text")
            valid = False
    if currYml.get('outputs'):
        if type(currYml['outputs']) is not dict:
            warnings.warn(f"{key} outputs is not an object")
            valid = False
        elif not all(d.get('purpose') and type(d['purpose']) is str
                     for d in currYml['outputs'].values()):
            warnings.warn(f"{key} output object is missing purpose text")
            valid = False
    if currYml.get('nameValue'):
        if type(currYml['nameValue']) is not dict:
            warnings.warn(f"{key} nameValue is not an object")
            valid = False
        elif not all(d.get('purpose') and type(d['purpose']) is str
                     for d in currYml['nameValue'].values()):
            warnings.warn(f"{key} nameValue object is missing purpose text")
            valid = False
    if currYml.get('examples'):
        if type(currYml['examples']) is not list or \
                not all(type(d) is dict for d in currYml['examples']):
            warnings.warn(f"{key} examples is not a list of objects")
            valid = False
        elif not all(d.get('head') and type(d['head']) is str and
                     d.get('body') and type(d['body']) is str
                     for d in currYml['examples']):
            warnings.warn(f"{key} examples object is missing head or body text")
            valid = False
    if currYml.get('moreAbout'):
        if type(currYml['moreAbout']) is not list or \
                not all(type(d) is dict for d in currYml['moreAbout']):
            warnings.warn(f"{key} moreAbout is not a list of objects")
            valid = False
        elif not all(d.get('head') and type(d['head']) is str and
                     d.get('body') and type(d['body']) is str
                     for d in currYml['moreAbout']):
            warnings.warn(f"{key} moreAbout object is missing head or body text")
            valid = False
    if not currYml.get('name'):
        warnings.warn(f"{key} lacks field name. Filling with name={key}")
        currYml['name'] = key
    if currYml.get('type') == 'method':
        cls = currKey.split('.')
        if len(cls) != 2:
            warnings.warn(f"{key} is a method, but {currKey} is not well structured")
            valid = False
        elif not currYml.get('class'):
            cls = cls[0]
            cls = '.'.join([*currPkg, cls])
            warnings.warn(f"{key} lacks field class. Filling with class={cls}")
            currYml['class'] = cls
    return valid


def checkClass():
    global currYml, currPkg, currKey
    valid = True
    key = '.'.join([*currPkg, currKey])
    if currYml.get('properties'):
        if type(currYml['properties']) is not dict:
            warnings.warn(f"{key} properties is not a dict")
            valid = False
        elif not all(d.get('purpose') and type(d['purpose']) is str
                     for d in currYml['properties'].values()):
            warnings.warn(f"{key} property object is missing purpose text")
            valid = False
    if currYml.get('constructors'):
        if type(currYml['constructors']) is not dict:
            warnings.warn(f"{key} constructors is a dict")
            valid = False
        elif type(currYml['constructors'].get('syntax')) is not list or \
                not all(stx and type(stx) is str for stx in currYml['constructors']['syntax']):
            warnings.warn(f"{key} input object is missing syntax field")
            valid = False
    if currYml.get('methods'):
        if type(currYml['methods']) is not dict:
            warnings.warn(f"{key} methods is not a dict")
            valid = False
    if not currYml.get('name'):
        warnings.warn(f"{key} lacks field name. Filling with name={key}")
        currYml['name'] = key
    return valid


def parseTemplate(filePath, signature):
    global env, currYml, currIndex, currPkg, currKey, currCls
    # if not ymlData.get('type'):
    #     raise NameError(f"No type specified for {key}")
    if currYml.get('type') == 'class':
        template = env.get_template("class.m.j2")
        if not checkClass():
            warnings.warn(f"{currKey} is not a proper class, skipping.")
            return signature
        currIndex['classes'] += currKey
        if currCls and currCls != currKey:
            raise ValueError(
                f"Improper format. Class={currKey} found in @{currCls} class folder. Cannot properly parse method documentations")
        currCls = currKey
    elif currYml.get('type') in ('function', 'method'):
        template = env.get_template("function.m.j2")
        if not checkFunction():
            warnings.warn(f"{currKey} is not a proper function, skipping.")
            return signature
        if currYml['type'] == 'function':
            currIndex['functions'] += currKey
        # TODO: include methods to classes
    elif currYml.get('type') in ('package', 'index'):
        template = env.get_template("index.m.j2")
    else:
        warnings.warn(f"No type specified for {currKey}. Defaulting to function template.")
        template = env.get_template("function.m.j2")
    os.makedirs(os.path.dirname(filePath), exist_ok=True)
    with open(filePath, 'w') as fl:
        fl.write(template.render(currYml))
    if currYml.get('type') in ('function', 'method'):
        signature = appendSignature(signature)
    return signature


# Documentation generator
def appendSignature(signature):
    global currYml
    if signature is None:
        signature = {}
    return signature


def generateMFile(yml: dict, path: str, file: str, signature: dict = None):
    global currKey, currYml, docroot, fullYml, currDocPath
    fullYml = yml
    if len(yml) > 1:
        # Allow defining multiple yaml objects in the same file
        for currKey, currYml in yml.items():
            key = currKey.split('.')
            if len(key) == 1:
                currDocPath = path
                signature = parseTemplate(os.path.join(docroot, path, key[0] + '.m'), signature)
            elif len(key) == 2:
                # Check if constructor
                currDocPath = os.path.join(path, key[0])
                if key[0] == key[1]:
                    signature = parseTemplate(os.path.join(docroot, path, key[0], 'constructor.m'), signature)
                else:
                    signature = parseTemplate(os.path.join(docroot, path, key[0], key[1] + '.m'), signature)
            else:
                warnings.warn(f"{currKey} must have either one index (function, class, package) or two indices ("
                              f"methods). Skipping")
                continue
    else:
        currKey, currYml = list(yml.items())[0]
        currDocPath = path
        # Should prioritize file for special cases like constructor
        signature = parseTemplate(os.path.join(docroot, path, file + '.m'), signature)
    return signature


def generateAllMFiles(path: str = "data"):
    global rootPath, currPath, currIndex, currPkg, currCls, htmlroot

    def getIndex(rPath: str = ''):
        nonlocal path
        ind = None
        if os.path.exists(file := os.path.join(path, rPath, 'index.yaml')) or \
                os.path.exists(file := os.path.join(path, rPath, 'index.yml')):
            with open(os.path.join(currPath, file)) as fl:
                ind = yaml.load(fl, yaml.Loader)
            if not ind or ind.get('type') not in ('index', 'package'):
                ind = None
        if not ind:
            if rPath:
                ind = {'type': 'package',
                       'packages': {},
                       'classes': {},
                       'functions': {}}
            else:
                ind = {'type': 'index',
                       'packages': {},
                       'classes': {},
                       'functions': {}}
        return ind

    def updateIndex(ind: dict, val: dict, p: list):
        if not p:
            return val
        p2 = p.pop(0)
        i = ind['packages'].get(p2)
        i = updateIndex(i, val, p)
        ind['packages'][p2] = i
        return ind

    def generateIndexM(ind: dict, toc: ET.Element, rpath: str = '', pkgs: list[str] = None):
        global currPkg
        if pkgs is None:
            pkgs = []
        key = pkgs[-1] if pkgs else 'index'
        currPkg = pkgs
        generateMFile({key: ind}, rpath, 'index')
        for pkg in ind['packages'].keys():
            tocItem = ET.Element(tag='tocitem')
            tocItem.attrib['target'] = '/'.join(['ref', *pkgs, pkg, 'index.html'])
            tocItem.text = '+' + pkg
            toc.append(tocItem)
            generateIndexM(ind['packages'][pkg], toc, os.path.join(relPath, pkg), [*pkgs, pkg])
        for cls in ind['classes']:
            tocItem = ET.Element(tag='tocitem')
            # tocItem.attrib['target'] = '/'.join(['ref', *pkgs, cls, 'index.html'])
            tocItem.attrib['target'] = '/'.join(['ref', *pkgs, cls, cls + '.html'])
            tocItem.text = '@' + cls
            toc.append(tocItem)
        # TODO: include functions and methods as well

    rootPath = path
    signature = None
    # Toolbox reference index
    index = {}
    # Generate all regular documentations and index packages, classes and functions
    for currPath, dirs, files in os.walk(path):
        relPath = os.path.relpath(currPath, path)
        currPkg = relPath.split(os.sep) if relPath else []
        # If a new package is being handled update the current index
        if all(pkg.startswith('+') for pkg in currPkg):
            currIndex = getIndex(relPath)
        # If it is a class folder, mark curCls
        currCls = currPkg[-1].removeprefix('@') if currPkg[-1].startswith('@') else ''
        # Extract the current/parent package
        currPkg = [pkg.removeprefix('+') for pkg in currPkg if pkg.startswith('+')]
        # Remove non-compliant folders (not +packages or @classes)
        for d in [d for d in dirs if not d.startswith(('+', '@'))]:
            # Have to use remove to make sure original dirs is changed
            dirs.remove(d)
        # Make sure @classes are iterated before +packages so that currIndex corresponds to the current package
        dirs.sort(reverse=True)
        # Iterate over documentation files
        # Only handle yaml files and ignore special case of index
        for file, fileName in [(file, fileName) for file in files
                               for fileName, ext in os.path.splitext(file)
                               if ext in ('.yaml', '.yml') and fileName != 'index']:
            with open(os.path.join(currPath, file)) as fl:
                yml = yaml.load(fl, yaml.Loader)
            signature = generateMFile(yml, relPath.replace('+', '').replace('@', ''), fileName, signature)

        # Update index with the currIndex updated by the current loop
        # In the case of @classes, this updates the parent package
        updateIndex(index, currIndex, currPkg)

    # Get current documentation toc
    tocFile = os.path.join(htmlroot, 'helptoc.xml')
    tree = ET.parse(tocFile)
    tocIndex = ET.Element(tag='tocitem')
    tocIndex.attrib['id'] = 'ref'
    tocIndex.attrib['target'] = 'ref/index.html'
    tocIndex.text = 'References'
    # Generate index files of the documentation and all packages
    generateIndexM(index, tocIndex)
    # Overwrite toc item
    if oldIndex := tree.find("./tocitem/tocitem[@id='ref']"):
        oldIndex.clear()
        # oldIndex.tag = tocIndex.tag
        oldIndex.attrib['id'] = tocIndex.attrib['id']
        oldIndex.attrib['target'] = tocIndex.attrib['target']
        oldIndex.text = tocIndex.text
        for tocItem in tocIndex.iterfind('./'):
            oldIndex.append(tocItem)
    else:
        tree.find("./tocitem").append(tocIndex)
    tree.write(tocFile)

    return signature


def main(root, path="data"):
    global env
    os.chdir(root)
    env.filters['asOutputs'] = asOutputs
    env.filters['asInputs'] = asInputs
    env.filters['getDataType'] = getDataType
    env.filters['getTypes'] = getTypes
    env.filters['getFQN'] = getFQN
    env.filters['getPurpose'] = getPurpose
    env.filters['isHTML'] = isHTML
    env.filters['getHRef'] = getHRef
    env.globals['getPackage'] = getPackage
    generateAllMFiles(path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='''
    Generate Matlab .m documentation files from the .j2 templates and .yaml data.\n
    Execute publish in Matlab with the StyleSheet.xsl (or simply publishAll.m) to generate the appropriate .html files.
    ''')
    parser.add_argument('root', default=os.getcwd(), nargs='?',
                        help='Root directory to calculate the relative paths from [default: %(default)s]')
    parser.add_argument('--doc', dest='doc', default=docroot,
                        help='Destination of the generated .m doc files [default: %(default)s]')
    parser.add_argument('--html', dest='html', default=htmlroot,
                        help='Location of the html root folder [default: %(default)s]')
    parser.add_argument('--data', dest='data', default='data',
                        help='Root directory containing the .yaml data files. Nested structure will be replicated in doc location [default: %(default)s]')
    try:
        args = parser.parse_args()
        docroot = args.doc
        htmlroot = args.html
        main(args.root, path=args.data)
    except SystemExit as e:
        if e.code == 0:
            pass
        else:
            raise
