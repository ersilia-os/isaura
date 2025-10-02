from isaura.manager.migrate import IsauraMover, IsauraCopy, IsauraRemover

IsauraCopy("eos3b5d", "v1", "gardp-project").copy()

# IsauraMover("eos3b5d", "v1", "gardp-project").move()

# IsauraRemover("eos3b5d", "v1", "gardp-project").remove()