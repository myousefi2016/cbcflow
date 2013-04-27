
def import_module_under_test(modulename):
    "Import module under test and print some info about it."
    mod = __import__(modulename)
    print("Running tests with %s version %s, date %s, imported from\n%s" % (
          modulename, mod.__version__, mod.__date__, mod.__file__))
    return mod

def init_test(name, modulename="headflow"):
    # Do nothing if the called module was imported
    if name == "__main__":
        # Use installed version from path if "i" is given as argument, otherwise import locally
        import sys
        if "i" in sys.argv[1:]:
            sys.argv.remove("i")
        else:
            sys.path.insert(0, "../../site-packages")

        mut = import_module_under_test(modulename)

