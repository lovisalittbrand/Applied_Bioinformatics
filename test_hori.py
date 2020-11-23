import argparse

def parseArguments():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument("creditMom", help="Credit mom.", type=str)
    parser.add_argument("creditDad", help="Credit dad.", type=float)
    parser.add_argument("debtMom", help="Debt mom.", type=float)

    # Optional arguments
    parser.add_argument("-dD", "--debtDad", help="Debt dad.", type=float, default=1000.)
    parser.add_argument("-s", "--salary", help="Debt dad.", type=float, default=2000.)
    parser.add_argument("-b", "--bonus", help="Debt dad.", type=float, default=0.)

    # Print version
    parser.add_argument("--version", action="version", version='%(prog)s - Version 1.0')

    # Parse arguments
    args = parser.parse_args()

    return args

def example(credit_mom, credit_dad, debt_mom, debt_dad = 1000, salary = 2000, bonus = 0):
    total_gain = salary + credit_dad + credit_mom + bonus
    total_loss = debt_dad + debt_mom

    return total_gain - total_loss


# Parse the arguments
args = parseArguments()

# Raw print arguments
print("You are running the script with arguments: ")
print(args.creditMom)
#for a in args.__dict__:
#    print(str(a))

# Run function
#print(example(args.creditMom, args.creditDad, args.debtMom, args.debtDad, args.salary, args.bonus))