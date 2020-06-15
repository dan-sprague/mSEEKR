if not args.dir.endswith('/'):
    args.dir+='/'
newDir = f'{args.dir}{args.qPrefix}_{args.nPrefix}/'
if not os.path.exists(newDir):
    os.mkdir(newDir)
else:
    flag = True
    while flag:
        usrIN = input(f'Directory {newDir} exists, continue? y/n: ').strip().lower()
        if usrIN == 'y':
            flag = False
        elif usrIN == 'n':
            print('Initiating self-destruct sequence')
            sys.exit()
        else:
            print('Please enter y or n')
