const restrictedSyntax = (...restrictions) => ({
  create(context) {
    return Object.fromEntries(
      restrictions.map(({ selector, message }) => [
        selector,
        node => {
          context.report({ node, message })
        },
      ]),
    )
  },
})

export default {
  meta: { name: 'local' },
  rules: {
    'no-spread-arguments': restrictedSyntax({
      selector: ':matches(CallExpression, NewExpression) > SpreadElement',
      message: 'Spread arguments overflow the stack on big arrays; use a loop',
    }),
  },
}
